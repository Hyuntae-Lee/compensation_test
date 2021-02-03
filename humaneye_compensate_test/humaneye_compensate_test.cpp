// humaneye_compensate_test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>
#include <vector>
#include <tuple>
#include <sstream>
#include <fstream>
#include <functional>
#include <opencv\cv.h>
#include <opencv2\imgproc.hpp>
#include <opencv2\imgcodecs.hpp>
#include <iostream>
#include <math.h>

#define MPI 3.141592
#define	DEGREE2RADIAN(_X) ((_X) / 180 * MPI)
#define RADIAN2DEGREE(_X) ((180.0/MPI) * _X)

using namespace std;

void meridian_fitting_test();
void fill_sca_data_hrk();
void fill_sca_data_cvs_002();
void fill_sca_data_cvs_004();
void fill_sca_data_cvs_005();
void fill_sample_data();
void calc_sca_ellipse(double CVS_meridian[3], double CVS_sca[4]);

void humaneye_compensate_test();
void humaneye_compensate(double sph, double cyl, double *tsph, double *tcyl);
auto getPowerForMeridians(double s, double c, double radian)->tuple<double, double, double>;

vector<pair<double, double>> s_sampleDataList;
vector<tuple<double, double, double>> s_scaDataList;

int main()
{
	vector<tuple<double, double, double>> meridianList = {
		{ -3.00, -2.45, -2.75 },
		{ -3.02, -2.42, -2.79 },
		{ -3.03, -2.43, -2.75 },
		{ -3.10, -2.55, -2.36 },
		{ -3.11, -2.53, -2.38 },
		{ -3.08, -2.56, -2.39 },
		{ -0.81, -1.02, -1.10 },
		{ -0.74, -0.94, -1.01 },
		{ -0.67, -0.89, -0.93 },
		{ -0.96, -0.82, -0.77 },
		{ -0.94, -0.88, -0.76 },
		{ -1.11, -0.96, -0.90 },
		{ -5.94, -6.19, -6.37 },
		{ -5.91, -6.17, -6.29 },
		{ -5.90, -6.14, -6.25 },
		{ -6.14, -7.16, -7.21 },
		{ -6.13, -7.13, -7.19 },
		{ -6.14, -7.14, -7.25 },
		{ -4.90, -5.77, -5.61 },
		{ -4.88, -5.75, -5.59 },
		{ -4.88, -5.75, -5.62 },
		{ -2.81, -3.62, -3.38 },
		{ -2.82, -3.63, -3.39 },
		{ -2.81, -3.60, -3.40 },
		{ -4.22, -3.89, -4.37 },
		{ -4.09, -3.76, -4.27 },
		{ -4.07, -3.74, -4.22 },
		{ -3.03, -3.48, -2.93 },
		{ -3.06, -3.47, -2.95 },
		{ -3.11, -3.50, -3.01 },
	};

	for (auto item : meridianList) {
		double CVS_meridian[3];
		
		CVS_meridian[0] = get<0>(item);
		CVS_meridian[1] = get<1>(item);
		CVS_meridian[2] = get<2>(item);
		
		double CVS_sca[4] = { 0.0, };

		calc_sca_ellipse(CVS_meridian, CVS_sca);

		auto s = CVS_sca[0];
		auto c = CVS_sca[1];
		auto a = (int)RADIAN2DEGREE(CVS_sca[2]);

		stringbuf strBuf;
		ostream os(&strBuf);

		os << "[";
		os << CVS_meridian[0];
		os << ", ";
		os << CVS_meridian[1];
		os << ", ";
		os << CVS_meridian[2];
		os << "] -> [";
		os << s;
		os << ", ";
		os << c;
		os << ", ";
		os << a;
		os << "]\r\n";

		cout << strBuf.str();
	}

	return 0;
}

void meridian_fitting_test()
{
	//fill_sca_data_hrk();
	//fill_sca_data_cvs_002();
	//fill_sca_data_cvs_004();
	fill_sca_data_cvs_005();

	stringbuf strBuf;
	ostream os(&strBuf);

	for (auto item : s_scaDataList) {
		auto s = get<0>(item);
		auto c = get<1>(item);
		auto a = get<2>(item);
		auto radian = DEGREE2RADIAN(a);

		tuple<double, double, double> meridians;

		if (s >= 0.0 || c >= 0.0) {
			meridians = make_tuple(0.0, 0.0, 0.0);
		}
		else {
			meridians = getPowerForMeridians(s, c, radian);
		}

		// fix the valid digits
		const int kStrBufSize = 50;

		char cMer000[kStrBufSize] = { 0, };
		char cMer060[kStrBufSize] = { 0, };
		char cMer120[kStrBufSize] = { 0, };

		sprintf_s(cMer000, kStrBufSize, "-%02.03f", get<0>(meridians));
		sprintf_s(cMer060, kStrBufSize, "-%02.03f", get<1>(meridians));
		sprintf_s(cMer120, kStrBufSize, "-%02.03f", get<2>(meridians));

		os << "[";
		os << s;
		os << ", ";
		os << c;
		os << ", ";
		os << a;
		os << "]\t=>\t[";
		os << cMer000;
		os << ", ";
		os << cMer060;
		os << ", ";
		os << cMer120;
		os << "]\n";
	}

	string contents = strBuf.str();

	ofstream ofs("D:/TTT/cvs_result.txt");
	ofs << contents;
	ofs.close();
}

auto getPowerForMeridians(double s, double c, double radian)->tuple<double, double, double>
{
	auto theta = radian;
	auto a = s;
	auto b = (s + c);

	auto thetaForAxis = 0.0;
	// degree 0
	thetaForAxis = -theta;
	auto powerFor000 = a * b * sqrt(1 / (pow(b, 2) * pow(cos(thetaForAxis), 2) + pow(a, 2) * pow(sin(thetaForAxis), 2)));
	// degree 60
	thetaForAxis = MPI / 3.0 - theta;
	auto powerFor060 = a * b * sqrt(1 / (pow(b, 2) * pow(cos(thetaForAxis), 2) + pow(a, 2) * pow(sin(thetaForAxis), 2)));
	// degree 120
	thetaForAxis = MPI * 2.0 / 3.0 - theta;
	auto powerFor120 = a * b * sqrt(1 / (pow(b, 2) * pow(cos(thetaForAxis), 2) + pow(a, 2) * pow(sin(thetaForAxis), 2)));

	return make_tuple(powerFor000, powerFor060, powerFor120);
}

void calc_sca_ellipse(double CVS_meridian[3], double CVS_sca[4])
{
	vector<cv::Point2f> ptList;

	float x_for_000 = (float)(-CVS_meridian[0]);
	float y_for_000 = 0.0f;
	float x_for_060 = (float)(-CVS_meridian[1] * cos(DEGREE2RADIAN(60.0)));
	float y_for_060 = (float)(-CVS_meridian[1] * sin(DEGREE2RADIAN(60.0)));
	float x_for_120 = (float)(-CVS_meridian[2] * cos(DEGREE2RADIAN(120.0)));
	float y_for_120 = (float)(-CVS_meridian[2] * sin(DEGREE2RADIAN(120.0)));
	float x_for_180 = -x_for_000;
	float y_for_180 = -y_for_000;
	float x_for_240 = -x_for_060;
	float y_for_240 = -y_for_060;
	float x_for_300 = -x_for_120;
	float y_for_300 = -y_for_120;

	ptList.push_back(cv::Point2f(x_for_000, y_for_000));
	ptList.push_back(cv::Point2f(x_for_060, y_for_060));
	ptList.push_back(cv::Point2f(x_for_120, y_for_120));
	ptList.push_back(cv::Point2f(x_for_180, y_for_180));
	ptList.push_back(cv::Point2f(x_for_240, y_for_240));
	ptList.push_back(cv::Point2f(x_for_300, y_for_300));

	cv::RotatedRect retRect;
	try {
		retRect = cv::fitEllipse(ptList);
	}
	catch (cv::Exception e) {
		auto strmsg = e.msg;
	}

	auto sph1 = retRect.size.width / 2.0f;
	auto sph2 = retRect.size.height / 2.0f;
	auto cyl = sph1 - sph2;

	CVS_sca[0] = sph1;
	CVS_sca[1] = sph2 - sph1;
	CVS_sca[2] = DEGREE2RADIAN(retRect.angle);
}

void fill_sca_data_hrk()
{
	s_scaDataList = {
		make_tuple(-0.74, -0.11, 73),
		make_tuple(-0.74, -0.08, 91),
		make_tuple(-0.75, -0.08, 91),
		make_tuple(-0.98, -0.17, 25),
		make_tuple(-0.96, -0.09, 17),
		make_tuple(-0.95, -0.12, 21),
		make_tuple(-0.69, -0.85, 81),
		make_tuple(-0.72, -0.84, 80),
		make_tuple(-0.63, -0.86, 80),
		make_tuple(-1.10, -0.32, 84),
		make_tuple(-1.05, -0.33, 88),
		make_tuple(-0.96, -0.36, 87),
		make_tuple(-1.80, -2.18, 174),
		make_tuple(-1.82, -2.16, 174),
		make_tuple(-1.81, -2.12, 175),
		make_tuple(-1.53, -1.68, 173),
		make_tuple(-1.66, -1.67, 173),
		make_tuple(-1.64, -1.68, 174),
		make_tuple(-0.23, -0.22, 112),
		make_tuple(-0.28, -0.20, 114),
		make_tuple(-0.12, -0.23, 111),
		make_tuple(-0.09, -0.34, 65),
		make_tuple(-0.09, -0.36, 67),
		make_tuple(-0.02, -0.35, 72),
		make_tuple(-5.00, -0.68, 104),
		make_tuple(-4.97, -0.68, 104),
		make_tuple(-4.98, -0.68, 104),
		make_tuple(-5.65, -0.22, 62),
		make_tuple(-5.59, -0.24, 68),
		make_tuple(-5.59, -0.24, 73),
		make_tuple(-0.99, -0.37, 5),
		make_tuple(-0.04, -0.36, 180),
		make_tuple(-0.08, -0.34, 1),
		make_tuple(-0.01, -0.40, 7),
		make_tuple(-0.14, -0.46, 6),
		make_tuple(-0.14, -0.36, 180),
		make_tuple(-0.82, -0.33, 8),
		make_tuple(-0.85, -0.37, 8),
		make_tuple(-0.87, -0.36, 8),
		make_tuple(-0.57, -0.43, 40),
		make_tuple(-0.58, -0.44, 30),
		make_tuple(-0.50, -0.38, 44),
		make_tuple(-0.59, -0.42, 179),
		make_tuple(-0.56, -0.44, 179),
		make_tuple(-0.52, -0.46, 180),
		make_tuple(-0.77, -0.09, 75),
		make_tuple(-0.82, -0.11, 70),
		make_tuple(-0.80, -0.18, 75),
		make_tuple(-0.25, -0.79, 82),
		make_tuple(-0.20, -0.81, 82),
		make_tuple(-0.20, -0.82, 81),
		make_tuple(0.14, -0.83, 93),
		make_tuple(0.08, -0.85, 89),
		make_tuple(0.13, -0.86, 89),
		make_tuple(-1.72, -2.39, 163),
		make_tuple(-1.81, -2.17, 162),
		make_tuple(-1.76, -2.22, 162),
		make_tuple(-5.14, -0.42, 171),
		make_tuple(-5.12, -0.41, 172),
		make_tuple(-5.15, -0.41, 172),
		make_tuple(-4.71, -1.22, 177),
		make_tuple(-4.69, -1.20, 176),
		make_tuple(-4.67, -1.22, 176),
		make_tuple(-2.69, -1.28, 170),
		make_tuple(-2.68, -1.29, 169),
		make_tuple(-2.65, -1.29, 170),
		make_tuple(-2.44, -0.76, 78),
		make_tuple(-2.47, -0.72, 77),
		make_tuple(-2.43, -0.73, 76),
		make_tuple(-2.12, -0.90, 101),
		make_tuple(-2.10, -0.93, 101),
		make_tuple(-2.11, -0.91, 101),
		make_tuple(-1.30, -0.69, 88),
		make_tuple(-1.29, -0.65, 89),
		make_tuple(-1.31, -0.59, 89),
		make_tuple(-0.73, -1.56, 81),
		make_tuple(-0.90, -1.61, 81),
		make_tuple(-0.92, -1.60, 82),
		make_tuple(-6.65, -0.49, 118),
		make_tuple(-6.78, -0.42, 125),
		make_tuple(-6.52, -0.52, 106),
		make_tuple(-5.63, -0.85, 12),
		make_tuple(-5.59, -0.92, 10),
		make_tuple(-5.57, -0.95, 9),
		make_tuple(-0.46, -1.44, 112),
		make_tuple(-0.50, -1.26, 104),
		make_tuple(-0.29, -1.32, 104),
		make_tuple(-0.96, -1.03, 67),
		make_tuple(-0.97, -1.03, 67),
		make_tuple(-0.92, -1.06, 65),
		make_tuple(-1.66, -0.73, 144),
		make_tuple(-1.91, -0.77, 147),
		make_tuple(-1.68, -0.77, 145),
		make_tuple(-1.56, -0.45, 68),
		make_tuple(-1.50, -0.41, 66),
		make_tuple(-1.47, -0.38, 68),
		make_tuple(-5.38, -0.27, 79),
		make_tuple(-5.40, -0.28, 81),
		make_tuple(-5.38, -0.29, 82),
		make_tuple(-5.85, -0.19, 2),
		make_tuple(-5.75, -0.15, 4),
		make_tuple(-5.79, -0.12, 10),
		make_tuple(-3.66, -1.10, 168),
		make_tuple(-3.64, -1.09, 167),
		make_tuple(-3.63, -1.10, 168),
		make_tuple(-2.99, -1.76, 7),
		make_tuple(-3.00, -1.78, 7),
		make_tuple(-2.98, -1.80, 7),
		make_tuple(-6.45, -1.62, 177),
		make_tuple(-6.45, -1.68, 177),
		make_tuple(-6.43, -1.69, 176),
		make_tuple(-5.12, -1.04, 5),
		make_tuple(-5.07, -1.15, 4),
		make_tuple(-5.10, -1.17, 5),
		make_tuple(-3.86, -0.52, 56),
		make_tuple(-3.84, -0.52, 53),
		make_tuple(-3.85, -0.44, 47),
		make_tuple(-3.76, -0.31, 104),
		make_tuple(-3.76, -0.28, 107),
		make_tuple(-3.71, -0.35, 114),
		make_tuple(-5.99, -0.46, 6),
		make_tuple(-5.97, -0.50, 7),
		make_tuple(-5.97, -0.54, 7),
		make_tuple(-6.29, -1.71, 180),
		make_tuple(-6.29, -1.84, 1),
		make_tuple(-6.27, -1.90, 180),
		make_tuple(-3.27, -0.54, 72),
		make_tuple(-3.23, -0.54, 75),
		make_tuple(-3.20, -0.56, 70),
		make_tuple(-3.28, -0.46, 106),
		make_tuple(-3.31, -0.42, 109),
		make_tuple(-3.28, -0.45, 109),
		make_tuple(-0.89, -0.47, 95),
		make_tuple(-0.83, -0.48, 94),
		make_tuple(-0.90, -0.47, 97),
		make_tuple(-0.36, -0.91, 85),
		make_tuple(-0.34, -0.90, 85),
		make_tuple(-0.33, -0.97, 85),
		make_tuple(-0.23, -0.54, 101),
		make_tuple(-0.18, -0.62, 94),
		make_tuple(-0.20, -0.59, 98),
		make_tuple(-0.34, -0.20, 160),
		make_tuple(-0.34, -0.21, 160),
		make_tuple(-0.30, -0.21, 155),
		make_tuple(-3.41, -0.49, 52),
		make_tuple(-3.41, -0.49, 49),
		make_tuple(-3.43, -0.46, 48),
		make_tuple(-2.93, -0.56, 148),
		make_tuple(-2.97, -0.54, 147),
		make_tuple(-2.97, -0.56, 149),
		make_tuple(-5.19, -0.53, 46),
		make_tuple(-5.17, -0.54, 47),
		make_tuple(-5.18, -0.50, 49),
		make_tuple(-4.03, -1.38, 169),
		make_tuple(-4.03, -1.36, 168),
		make_tuple(-3.99, -1.37, 169),
		make_tuple(-0.18, -1.20, 161),
		make_tuple(-0.18, -1.18, 161),
		make_tuple(-0.22, -1.14, 161),
		make_tuple(-0.03, -1.96, 7),
		make_tuple(-0.16, -1.90, 7),
		make_tuple(-0.26, -1.86, 8),
		make_tuple(-2.16, -0.76, 26),
		make_tuple(-2.19, -0.76, 27),
		make_tuple(-2.15, -0.77, 26),
		make_tuple(-2.77, -0.33, 136),
		make_tuple(-2.78, -0.31, 133),
		make_tuple(-2.77, -0.31, 137),
		make_tuple(-5.81, -1.23, 172),
		make_tuple(-5.75, -1.30, 173),
		make_tuple(-5.74, -1.35, 172),
		make_tuple(-6.61, -0.68, 173),
		make_tuple(-6.57, -0.73, 176),
		make_tuple(-6.63, -0.75, 177),
	};
}

void fill_sca_data_cvs_002()
{
	s_scaDataList = {
		make_tuple(-0.24, -0.32, 114),
		make_tuple(-0.33, -0.38, 109),
		make_tuple(-0.26, -0.33, 107),
		make_tuple(-0.39, -0.40, 11),
		make_tuple(-0.50, -0.19, 21),
		make_tuple(-0.49, -0.23, 16),
		make_tuple(-0.09, -0.44, 89),
		make_tuple(0.11, -0.53, 106),
		make_tuple(-0.11, -0.36, 97),
		make_tuple(0.04, -0.53, 59),
		make_tuple(0.28, -0.70, 53),
		make_tuple(0.05, -0.51, 56),
		make_tuple(-0.93, -2.77, 10),
		make_tuple(-1.08, -2.77, 1),
		make_tuple(-0.74, -2.77, 10),
		make_tuple(-0.99, -1.76, 179),
		make_tuple(-0.88, -1.80, 176),
		make_tuple(-0.59, -1.87, 173),
		make_tuple(0.12, -0.53, 130),
		make_tuple(0.50, -0.87, 126),
		make_tuple(0.19, -0.37, 138),
		make_tuple(0.28, -0.24, 61),
		make_tuple(0.53, -0.39, 66),
		make_tuple(0.55, -0.35, 41),
		make_tuple(-4.30, -0.83, 106),
		make_tuple(-4.62, -0.66, 101),
		make_tuple(-4.41, -1.06, 82),
		make_tuple(-4.50, -0.67, 50),
		make_tuple(-5.20, -0.14, 66),
		make_tuple(-5.11, -0.18, 52),
		make_tuple(-0.37, -0.12, 160),
		make_tuple(-0.34, -0.25, 143),
		make_tuple(-0.51, -0.05, 113),
		make_tuple(-0.46, -0.34, 15),
		make_tuple(-0.52, -0.29, 0),
		make_tuple(-0.58, -0.25, 178),
		make_tuple(-0.63, -0.45, 149),
		make_tuple(-0.78, -0.28, 157),
		make_tuple(-0.68, -0.49, 152),
		make_tuple(-0.56, -0.18, 25),
		make_tuple(-0.63, -0.07, 14),
		make_tuple(-0.59, -0.16, 38),
		make_tuple(-0.57, -0.44, 18),
		make_tuple(-0.38, -0.42, 16),
		make_tuple(-0.46, -0.41, 21),
		make_tuple(-0.78, -0.16, 8),
		make_tuple(-0.65, -0.11, 24),
		make_tuple(-0.57, -0.20, 47),
		make_tuple(0.14, -0.78, 100),
		make_tuple(-0.10, -0.40, 85),
		make_tuple(-0.05, -0.59, 93),
		make_tuple(0.24, -0.47, 79),
		make_tuple(0.27, -0.46, 81),
		make_tuple(0.29, -0.63, 77),
		make_tuple(-1.07, -1.83, 0),
		make_tuple(-1.17, -1.97, 7),
		make_tuple(-1.13, -2.07, 10),
		make_tuple(-4.37, -0.59, 148),
		make_tuple(-4.46, -0.85, 148),
		make_tuple(-4.47, -0.67, 159),
		make_tuple(-4.22, -2.60, 15),
		make_tuple(-2.56, -2.77, 11),
		make_tuple(-2.16, -2.77, 9),
		make_tuple(-1.55, -1.33, 166),
		make_tuple(-1.21, -2.77, 161),
		make_tuple(-0.84, -2.70, 164),
		make_tuple(-2.00, -0.82, 66),
		make_tuple(-1.96, -0.89, 61),
		make_tuple(-1.90, -0.96, 69),
		make_tuple(-1.44, -1.22, 99),
		make_tuple(-1.52, -0.88, 104),
		make_tuple(-1.51, -1.23, 101),
		make_tuple(-1.08, -0.73, 57),
		make_tuple(-0.22, -0.94, 92),
		make_tuple(-0.50, -0.52, 98),
		make_tuple(-1.07, -0.73, 84),
		make_tuple(-0.36, -1.06, 69),
		make_tuple(-0.81, -0.73, 74),
		make_tuple(-6.04, -0.29, 45),
		make_tuple(-6.08, -0.75, 32),
		make_tuple(-6.05, -0.46, 33),
		make_tuple(-4.57, -0.78, 20),
		make_tuple(-4.70, -0.99, 14),
		make_tuple(-4.72, -1.20, 6),
		make_tuple(0.06, -0.93, 123),
		make_tuple(0.24, -1.10, 111),
		make_tuple(0.08, -0.91, 116),
		make_tuple(-0.47, -0.89, 68),
		make_tuple(-0.49, -0.88, 45),
		make_tuple(-0.40, -0.65, 57),
		make_tuple(-0.28, -0.58, 133),
		make_tuple(-0.28, -0.60, 126),
		make_tuple(-0.26, -0.56, 121),
		make_tuple(-0.45, -0.46, 57),
		make_tuple(-0.72, -0.20, 47),
		make_tuple(-0.62, -0.39, 67),
		make_tuple(-5.84, -0.25, 103),
		make_tuple(-5.12, -0.74, 94),
		make_tuple(-4.98, -1.03, 107),
		make_tuple(-5.66, -0.45, 21),
		make_tuple(-5.52, -0.25, 43),
		make_tuple(-5.38, -0.45, 67),
		make_tuple(-3.45, -1.21, 164),
		make_tuple(-3.75, -0.98, 170),
		make_tuple(-3.49, -1.29, 167),
		make_tuple(-1.99, -2.67, 3),
		make_tuple(-2.09, -2.60, 3),
		make_tuple(-2.08, -2.70, 5),
		make_tuple(-5.44, -1.54, 8),
		make_tuple(-5.44, -1.55, 6),
		make_tuple(-5.68, -1.31, 5),
		make_tuple(-4.47, -1.21, 165),
		make_tuple(-4.54, -1.12, 165),
		make_tuple(-4.82, -1.01, 173),
		make_tuple(-4.07, -0.68, 31),
		make_tuple(-3.98, -1.17, 34),
		make_tuple(-4.07, -0.60, 30),
		make_tuple(-3.52, -1.19, 118),
		make_tuple(-3.53, -1.22, 125),
		make_tuple(-3.57, -1.19, 117),
		make_tuple(-5.28, -1.71, 36),
		make_tuple(-5.13, -1.86, 57),
		make_tuple(-6.13, -0.86, 34),
		make_tuple(-6.21, -0.78, 180),
		make_tuple(-5.24, -1.75, 177),
		make_tuple(-6.77, -0.22, 1),
		make_tuple(-2.60, -1.48, 61),
		make_tuple(-2.51, -1.59, 53),
		make_tuple(-2.43, -1.50, 58),
		make_tuple(-2.15, -1.40, 134),
		make_tuple(-2.11, -1.44, 130),
		make_tuple(-1.91, -1.11, 124),
		make_tuple(0.11, -0.62, 104),
		make_tuple(-0.10, -0.54, 97),
		make_tuple(0.08, -0.70, 97),
		make_tuple(-0.01, -0.88, 64),
		make_tuple(-0.20, -0.72, 70),
		make_tuple(-0.05, -1.00, 72),
		make_tuple(-0.13, -0.51, 110),
		make_tuple(-0.15, -0.37, 100),
		make_tuple(-0.19, -0.38, 111),
		make_tuple(-0.24, -0.05, 54),
		make_tuple(-0.11, -0.24, 171),
		make_tuple(-0.09, -0.25, 176),
		make_tuple(-2.57, -1.76, 47),
		make_tuple(-1.73, -2.14, 37),
		make_tuple(-1.79, -2.28, 36),
		make_tuple(-1.16, -1.86, 148),
		make_tuple(-1.34, -1.12, 141),
		make_tuple(-1.28, -0.90, 148),
		make_tuple(-4.59, -1.87, 40),
		make_tuple(-4.53, -1.96, 25),
		make_tuple(-5.40, -1.06, 38),
		make_tuple(-3.26, -1.98, 165),
		make_tuple(-2.57, -2.77, 167),
		make_tuple(-2.66, -2.77, 172),
		make_tuple(-0.24, -0.79, 157),
		make_tuple(-0.37, -0.71, 161),
		make_tuple(-0.36, -0.74, 171),
		make_tuple(-0.03, -1.46, 176),
		make_tuple(-0.09, -1.36, 178),
		make_tuple(0.01, -1.56, 177),
		make_tuple(-1.51, -0.87, 33),
		make_tuple(-1.27, -1.52, 32),
		make_tuple(-1.36, -1.40, 34),
		make_tuple(-1.68, -0.64, 147),
		make_tuple(-1.39, -0.92, 154),
		make_tuple(-1.63, -0.97, 144),
		make_tuple(-5.44, -1.50, 175),
		make_tuple(-5.62, -1.12, 162),
		make_tuple(-5.58, -0.96, 171),
		make_tuple(-6.44, -0.54, 174),
		make_tuple(-6.39, -0.60, 175),
		make_tuple(-6.32, -0.67, 159),
	};
}

void fill_sca_data_cvs_004()
{
	s_scaDataList = {
		make_tuple(-0.47, -0.28, 97),
		make_tuple(-0.48, -0.02, 47),
		make_tuple(-0.50, -0.02, 115),
		make_tuple(-0.55, -0.19, 15),
		make_tuple(-0.52, -0.01, 34),
		make_tuple(-0.47, -0.06, 172),
		make_tuple(-0.11, -0.63, 120),
		make_tuple(-0.27, -0.45, 113),
		make_tuple(0.22, -0.03, 59),
		make_tuple(0.22, -0.81, 47),
		make_tuple(0.03, -0.70, 52),
		make_tuple(0.23, -0.02, 36),
		make_tuple(-0.52, -2.41, 23),
		make_tuple(-0.85, -1.24, 5),
		make_tuple(-0.69, -2.46, 165),
		make_tuple(-0.72, -1.00, 171),
		make_tuple(-0.80, -0.76, 174),
		make_tuple(-0.79, -0.79, 175),
		make_tuple(0.04, -0.64, 119),
		make_tuple(-0.31, -1.03, 119),
		make_tuple(0.20, -0.02, 57),
		make_tuple(0.28, -0.51, 71),
		make_tuple(0.27, -0.45, 71),
		make_tuple(0.24, -0.06, 49),
		make_tuple(-4.52, -0.66, 88),
		make_tuple(-4.27, -0.37, 64),
		make_tuple(-2.09, -1.35, 3),
		make_tuple(-4.74, -0.62, 172),
		make_tuple(-4.27, -0.81, 149),
		make_tuple(-2.12, -1.35, 164),
		make_tuple(0.20, -0.03, 58),
		make_tuple(-0.34, -0.20, 75),
		make_tuple(-0.56, -0.04, 99),
		make_tuple(-0.47, -0.02, 2),
		make_tuple(-0.63, -0.11, 156),
		make_tuple(-0.72, -0.23, 1),
		make_tuple(-0.82, -0.35, 152),
		make_tuple(0.38, -0.09, 104),
		make_tuple(-0.51, -0.02, 149),
		make_tuple(-0.58, -0.49, 26),
		make_tuple(0.36, -0.06, 60),
		make_tuple(-0.04, -0.07, 43),
		make_tuple(-0.49, -0.46, 156),
		make_tuple(-0.69, -0.14, 2),
		make_tuple(-0.73, -0.22, 3),
		make_tuple(-0.61, -0.29, 38),
		make_tuple(-0.83, -0.07, 4),
		make_tuple(-0.81, -0.10, 31),
		make_tuple(0.15, -0.98, 76),
		make_tuple(-0.18, -0.64, 88),
		make_tuple(0.20, -0.04, 68),
		make_tuple(0.03, -0.59, 84),
		make_tuple(-0.11, -0.30, 83),
		make_tuple(0.22, -0.06, 42),
		make_tuple(0.36, -0.05, 117),
		make_tuple(-0.67, -2.77, 154),
		make_tuple(-0.92, -2.26, 160),
		make_tuple(0.40, -2.77, 124),
		make_tuple(-4.14, -0.81, 6),
		make_tuple(-4.19, -0.69, 174),
		make_tuple(0.40, -2.77, 149),
		make_tuple(-4.21, -0.34, 2),
		make_tuple(-4.91, -1.06, 168),
		make_tuple(-0.56, -0.07, 9),
		make_tuple(-1.46, -0.04, 178),
		make_tuple(-1.84, -2.64, 179),
		make_tuple(-1.60, -0.99, 83),
		make_tuple(-0.94, -0.06, 75),
		make_tuple(-0.49, -0.11, 59),
		make_tuple(-1.35, -0.72, 95),
		make_tuple(-0.55, -0.05, 99),
		make_tuple(-0.02, -0.09, 94),
		make_tuple(-0.47, -0.08, 117),
		make_tuple(-0.82, -0.67, 79),
		make_tuple(-0.52, -0.08, 67),
		make_tuple(-0.48, -0.06, 97),
		make_tuple(-1.04, -0.45, 74),
		make_tuple(-0.53, -0.04, 62),
		make_tuple(-5.46, -0.70, 43),
		make_tuple(-3.29, -0.12, 73),
		make_tuple(-3.67, -0.55, 51),
		make_tuple(-4.40, -0.79, 2),
		make_tuple(-2.41, -0.03, 1),
		make_tuple(-2.35, -0.09, 37),
		make_tuple(-0.13, -0.84, 109),
		make_tuple(-0.11, -0.71, 112),
		make_tuple(-0.07, -0.06, 101),
		make_tuple(-0.56, -0.66, 65),
		make_tuple(-0.03, -1.41, 61),
		make_tuple(-0.50, -0.02, 68),
		make_tuple(-0.24, -0.63, 119),
		make_tuple(0.21, -0.07, 121),
		make_tuple(-0.35, -0.62, 124),
		make_tuple(-0.50, -0.48, 59),
		make_tuple(0.39, -0.10, 57),
		make_tuple(-0.70, -0.27, 50),
		make_tuple(-4.30, -0.06, 58),
		make_tuple(-3.26, -0.13, 88),
		make_tuple(-3.31, -0.13, 61),
		make_tuple(-4.33, -0.11, 131),
		make_tuple(-3.31, -0.13, 146),
		make_tuple(-3.29, -0.20, 136),
		make_tuple(-1.47, -0.09, 60),
		make_tuple(-3.41, -1.27, 178),
		make_tuple(-3.50, -1.67, 171),
		make_tuple(-1.63, -1.69, 174),
		make_tuple(-1.89, -2.77, 173),
		make_tuple(-3.22, -1.92, 170),
		make_tuple(-5.85, -1.14, 167),
		make_tuple(-5.62, -1.37, 5),
		make_tuple(-4.28, -0.28, 169),
		make_tuple(-4.71, -1.39, 177),
		make_tuple(-4.46, -1.46, 171),
		make_tuple(-3.29, -0.15, 15),
		make_tuple(-1.58, -1.68, 177),
		make_tuple(0.40, -2.72, 153),
		make_tuple(-2.16, -1.19, 16),
		make_tuple(0.40, -2.47, 150),
		make_tuple(0.40, -2.57, 148),
		make_tuple(-1.66, -1.58, 117),
		make_tuple(-4.53, -2.02, 33),
		make_tuple(0.40, -2.42, 170),
		make_tuple(-4.30, -2.33, 35),
		make_tuple(-4.98, -2.01, 155),
		make_tuple(0.40, -2.59, 167),
		make_tuple(-4.97, -2.02, 155),
		make_tuple(-3.32, -0.24, 87),
		make_tuple(-1.53, -1.42, 75),
		make_tuple(-0.94, -0.04, 173),
		make_tuple(-1.98, -1.76, 128),
		make_tuple(-1.41, -0.05, 72),
		make_tuple(0.40, -2.35, 150),
		make_tuple(0.23, -0.09, 117),
		make_tuple(0.20, -0.05, 62),
		make_tuple(0.22, -0.07, 105),
		make_tuple(0.21, -0.05, 99),
		make_tuple(0.23, -1.89, 58),
		make_tuple(0.23, -0.08, 63),
		make_tuple(-0.16, -0.61, 93),
		make_tuple(-0.09, -0.05, 111),
		make_tuple(-0.16, -0.58, 103),
		make_tuple(-0.33, -0.14, 61),
		make_tuple(-0.49, -0.02, 38),
		make_tuple(-0.19, -0.17, 44),
		make_tuple(-1.90, -1.38, 160),
		make_tuple(-1.82, -1.38, 159),
		make_tuple(-1.87, -1.60, 177),
		make_tuple(-1.48, -0.11, 150),
		make_tuple(-1.47, -0.13, 150),
		make_tuple(-1.45, -0.16, 139),
		make_tuple(-4.19, -0.43, 52),
		make_tuple(-3.66, -0.29, 69),
		make_tuple(-3.30, -0.12, 104),
		make_tuple(-2.08, -1.59, 173),
		make_tuple(-1.65, -1.70, 137),
		make_tuple(-1.62, -1.67, 136),
		make_tuple(0.23, -1.90, 32),
		make_tuple(-0.39, -0.58, 156),
		make_tuple(-0.47, -0.52, 160),
		make_tuple(0.23, -1.91, 18),
		make_tuple(-0.10, -1.26, 177),
		make_tuple(-0.30, -1.09, 176),
		make_tuple(-1.24, -0.58, 31),
		make_tuple(-0.54, -0.05, 40),
		make_tuple(-1.16, -0.61, 28),
		make_tuple(-1.44, -0.25, 151),
		make_tuple(-0.56,  0.00, 93),
		make_tuple(-1.28, -0.45, 160),
		make_tuple(-5.31, -1.44, 169),
		make_tuple(-4.06, -0.53, 169),
		make_tuple(-5.25, -1.46, 173),
		make_tuple(-5.92, -1.07, 165),
		make_tuple(-4.55, -0.45, 161),
		make_tuple(-5.84, -1.15, 168),

	};
}

void fill_sca_data_cvs_005()
{
	s_scaDataList = {
		make_tuple(-0.40, -0.24, 83),
		make_tuple(-0.43, -0.14, 109),
		make_tuple(-0.50, -0.22, 93),
		make_tuple(-0.62, -0.08, 4),
		make_tuple(-0.58, -0.10, 151),
		make_tuple(-0.60, -0.04, 25),
		make_tuple(-0.39, -0.28, 74),
		make_tuple(-0.37, -0.32, 79),
		make_tuple(-0.38, -0.19, 69),
		make_tuple(-0.01, -0.67, 67),
		make_tuple(-0.02, -0.59, 77),
		make_tuple(0.27, -0.83, 69),
		make_tuple(-0.69, -2.01, 4),
		make_tuple(-0.59, -1.66, 178),
		make_tuple(-0.65, -1.64, 177),
		make_tuple(-0.74, -1.64, 172),
		make_tuple(-0.66, -1.29, 170),
		make_tuple(-0.70, -1.12, 179),
		make_tuple(0.01, -0.28, 158),
		make_tuple(-0.03, -0.40, 163),
		make_tuple(0.00, -0.10, 13),
		make_tuple(0.40, -0.62, 72),
		make_tuple(0.60, -0.47, 72),
		make_tuple(0.41, -0.63, 82),
		make_tuple(-4.60, -1.16, 81),
		make_tuple(-4.00, -1.07, 65),
		make_tuple(-4.34, -1.58, 109),
		make_tuple(-4.89, -0.79, 147),
		make_tuple(-5.19, -0.44, 166),
		make_tuple(-5.28, -0.45, 35),
		make_tuple(-0.59, -0.08, 152),
		make_tuple(-0.16, -0.35, 180),
		make_tuple(-0.43, -0.14, 150),
		make_tuple(-0.73, -0.09, 176),
		make_tuple(-0.29, -0.29, 79),
		make_tuple(-0.66, -0.06, 164),
		make_tuple(-0.94, -0.29, 164),
		make_tuple(-0.77, -0.22, 166),
		make_tuple(-0.73, -0.38, 160),
		make_tuple(-0.81, -0.25, 35),
		make_tuple(-0.51, -0.45, 57),
		make_tuple(-0.57, -0.34, 52),
		make_tuple(-0.60, -0.37, 3),
		make_tuple(-0.42, -0.26, 166),
		make_tuple(-0.70, -0.25, 7),
		make_tuple(-0.73, -0.21, 41),
		make_tuple(-0.64, -0.09, 85),
		make_tuple(-0.82, -0.08, 52),
		make_tuple(-0.18, -0.62, 83),
		make_tuple(-0.23, -0.35, 82),
		make_tuple(-0.23, -0.41, 81),
		make_tuple(0.40, -0.97, 83),
		make_tuple(0.38, -1.08, 83),
		make_tuple(0.27, -0.86, 83),
		make_tuple(-1.10, -1.25, 176),
		make_tuple(-0.95, -1.59, 15),
		make_tuple(-0.99, -1.22, 175),
		make_tuple(-3.74, -1.19, 173),
		make_tuple(-3.67, -1.78, 161),
		make_tuple(-2.97, -1.81, 165),
		make_tuple(-4.08, -0.43, 142),
		make_tuple(-1.84, -2.77, 16),
		make_tuple(-2.12, -2.77, 12),
		make_tuple(-1.49, -1.23, 3),
		make_tuple(-0.91, -2.77, 156),
		make_tuple(-1.13, -2.15, 170),
		make_tuple(-1.58, -1.38, 82),
		make_tuple(-1.80, -0.63, 69),
		make_tuple(-1.38, -0.91, 81),
		make_tuple(-1.54, -0.59, 106),
		make_tuple(-1.48, -0.63, 101),
		make_tuple(-1.31, -0.83, 103),
		make_tuple(-0.39, -1.62, 89),
		make_tuple(-0.57, -0.46, 86),
		make_tuple(-0.88, -0.50, 101),
		make_tuple(-0.96, -0.53, 71),
		make_tuple(-0.93, -0.40, 83),
		make_tuple(-1.14, -0.53, 71),
		make_tuple(-6.00, -0.17, 74),
		make_tuple(-5.57, -0.73, 30),
		make_tuple(-5.48, -0.79, 30),
		make_tuple(-4.29, -1.13, 14),
		make_tuple(-4.33, -1.48, 176),
		make_tuple(-4.45, -1.32, 173),
		make_tuple(-0.36, -0.31, 120),
		make_tuple(-0.23, -0.52, 125),
		make_tuple(-0.24, -0.53, 129),
		make_tuple(-0.17, -1.26, 63),
		make_tuple(-0.63, -0.55, 62),
		make_tuple(-0.70, -0.33, 53),
		make_tuple(0.08, -0.66, 84),
		make_tuple(-0.39, -0.51, 116),
		make_tuple(-0.29, -0.55, 114),
		make_tuple(-0.31, -0.72, 71),
		make_tuple(-0.73, -0.30, 61),
		make_tuple(-0.77, -0.26, 180),
		make_tuple(-4.55, -0.86, 61),
		make_tuple(-5.19, -0.16, 90),
		make_tuple(-4.76, -0.12, 174),
		make_tuple(-4.77, -1.04, 139),
		make_tuple(-5.06, -0.67, 155),
		make_tuple(-5.02, -1.26, 165),
		make_tuple(-3.12, -1.55, 175),
		make_tuple(-3.62, -1.22, 168),
		make_tuple(-3.55, -1.51, 164),
		make_tuple(-1.85, -2.67, 179),
		make_tuple(-1.66, -2.77, 3),
		make_tuple(-1.74, -2.77, 7),
		make_tuple(-5.24, -1.74, 18),
		make_tuple(-5.19, -1.80, 22),
		make_tuple(-4.89, -2.10, 15),
		make_tuple(-4.16, -2.30, 156),
		make_tuple(-4.13, -2.15, 151),
		make_tuple(-4.33, -1.87, 157),
		make_tuple(-3.42, -0.66, 70),
		make_tuple(-3.52, -0.60, 74),
		make_tuple(-3.47, -0.66, 68),
		make_tuple(-3.10, -1.25, 139),
		make_tuple(-3.05, -1.09, 134),
		make_tuple(-2.94, -1.36, 138),
		make_tuple(-6.01, -0.98, 52),
		make_tuple(-3.95, -1.54, 60),
		make_tuple(-6.30, -0.69, 44),
		make_tuple(-6.16, -0.83, 1),
		make_tuple(-2.72, -2.77, 175),
		make_tuple(-5.18, -1.81, 179),
		make_tuple(-1.83, -1.74, 53),
		make_tuple(-2.18, -1.38, 73),
		make_tuple(-2.10, -1.77, 70),
		make_tuple(-1.98, -1.76, 135),
		make_tuple(-1.94, -1.44, 148),
		make_tuple(-1.69, -1.72, 148),
		make_tuple(-0.26, -0.46, 105),
		make_tuple(-0.34, -0.37, 101),
		make_tuple(-0.16, -0.54, 106),
		make_tuple(-0.03, -0.95, 81),
		make_tuple(-0.07, -0.85, 79),
		make_tuple(-0.23, -0.76, 77),
		make_tuple(-0.35, -0.50, 114),
		make_tuple(-0.32, -0.43, 107),
		make_tuple(-0.29, -0.46, 114),
		make_tuple(-0.41, -0.07, 96),
		make_tuple(-0.36, -0.24, 113),
		make_tuple(-0.31, -0.17, 108),
		make_tuple(-2.14, -2.00, 65),
		make_tuple(-1.58, -1.82, 46),
		make_tuple(-1.59, -1.77, 45),
		make_tuple(-1.50, -0.54, 159),
		make_tuple(-1.02, -1.97, 150),
		make_tuple(-0.97, -2.04, 149),
		make_tuple(-3.71, -2.38, 54),
		make_tuple(-4.54, -1.08, 58),
		make_tuple(-3.95, -2.73, 49),
		make_tuple(-2.70, -2.41, 153),
		make_tuple(-1.88, -2.77, 163),
		make_tuple(-1.89, -2.77, 166),
		make_tuple(-0.44, -0.58, 166),
		make_tuple(-0.35, -0.79, 173),
		make_tuple(-0.34, -0.79, 166),
		make_tuple(-0.32, -1.05, 3),
		make_tuple(-0.34, -1.09, 2),
		make_tuple(-0.22, -1.23, 1),
		make_tuple(-1.43, -1.17, 35),
		make_tuple(-1.32, -0.95, 42),
		make_tuple(-1.19, -1.28, 38),
		make_tuple(-1.49, -0.67, 161),
		make_tuple(-1.70, -0.45, 146),
		make_tuple(-1.60, -0.57, 159),
		make_tuple(-6.00, -0.99, 163),
		make_tuple(-5.84, -1.15, 20),
		make_tuple(-5.93, -0.99, 13),
		make_tuple(-6.50, -0.49, 170),
		make_tuple(-5.60, -1.39, 155),
		make_tuple(-5.79, -1.20, 158),
	};
}

void humaneye_compensate_test()
{
	fill_sample_data();

	stringbuf strBuf;
	ostream os(&strBuf);

	for (auto item : s_sampleDataList) {
		double tsph, tcyl;
		humaneye_compensate(get<0>(item), get<1>(item), &tsph, &tcyl);

		os << tsph;
		os << ", ";
		os << tcyl;
		os << "\n";
	}

	string contents = strBuf.str();

	ofstream ofs("D:/TTT/cvs_result.txt");
	ofs << contents;
	ofs.close();
}

void humaneye_compensate(double sph, double cyl, double *tsph, double *tcyl)
{
	const double kCompRatio = 0.5;

	double sph1_raw = sph;
	double sph2_raw = sph + cyl;

	double sph1_on_graph = 0.0005 * pow(sph1_raw, 3) + 0.095 * pow(sph1_raw, 2) + 1.5601 * sph1_raw - 0.2966;
	double sph2_on_graph = 0.0005 * pow(sph2_raw, 3) + 0.095 * pow(sph2_raw, 2) + 1.5601 * sph2_raw - 0.2966;

	double sph1 = sph1_raw + (sph1_on_graph - sph1_raw) * kCompRatio;
	double sph2 = sph2_raw + (sph2_on_graph - sph2_raw) * kCompRatio;

	*tsph = sph1;
	*tcyl = sph2 - sph1;
}

void fill_sample_data()
{
	s_sampleDataList = {
		make_pair(-0.24, -0.32),
		make_pair(-0.33, -0.38),
		make_pair(-0.26, -0.33),
		make_pair(-0.39, -0.40),
		make_pair(-0.50, -0.19),
		make_pair(-0.49, -0.23),
		make_pair(-0.09, -0.44),
		make_pair(0.11, -0.53),
		make_pair(-0.11, -0.36),
		make_pair(0.04, -0.53),
		make_pair(0.28, -0.70),
		make_pair(0.05, -0.51),
		make_pair(-0.93, -2.77),
		make_pair(-1.08, -2.77),
		make_pair(-0.74, -2.77),
		make_pair(-0.99, -1.76),
		make_pair(-0.88, -1.80),
		make_pair(-0.59, -1.87),
		make_pair(0.12, -0.53),
		make_pair(0.50, -0.87),
		make_pair(0.19, -0.37),
		make_pair(0.28, -0.24),
		make_pair(0.53, -0.39),
		make_pair(0.55, -0.35),
		make_pair(-4.30, -0.83),
		make_pair(-4.62, -0.66),
		make_pair(-4.41, -1.06),
		make_pair(-4.50, -0.67),
		make_pair(-5.20, -0.14),
		make_pair(-5.11, -0.18),
		make_pair(-0.37, -0.12),
		make_pair(-0.34, -0.25),
		make_pair(-0.51, -0.05),
		make_pair(-0.46, -0.34),
		make_pair(-0.52, -0.29),
		make_pair(-0.58, -0.25),
		make_pair(-0.63, -0.45),
		make_pair(-0.78, -0.28),
		make_pair(-0.68, -0.49),
		make_pair(-0.56, -0.18),
		make_pair(-0.63, -0.07),
		make_pair(-0.59, -0.16),
		make_pair(-0.57, -0.44),
		make_pair(-0.38, -0.42),
		make_pair(-0.46, -0.41),
		make_pair(-0.78, -0.16),
		make_pair(-0.65, -0.11),
		make_pair(-0.57, -0.20),
		make_pair(0.14, -0.78),
		make_pair(-0.10, -0.40),
		make_pair(-0.05, -0.59),
		make_pair(0.24, -0.47),
		make_pair(0.27, -0.46),
		make_pair(0.29, -0.63),
		make_pair(-1.07, -1.83),
		make_pair(-1.17, -1.97),
		make_pair(-1.13, -2.07),
		make_pair(-4.37, -0.59),
		make_pair(-4.46, -0.85),
		make_pair(-4.47, -0.67),
		make_pair(-4.22, -2.60),
		make_pair(-2.56, -2.77),
		make_pair(-2.16, -2.77),
		make_pair(-1.55, -1.33),
		make_pair(-1.21, -2.77),
		make_pair(-0.84, -2.70),
		make_pair(-2.00, -0.82),
		make_pair(-1.96, -0.89),
		make_pair(-1.90, -0.96),
		make_pair(-1.44, -1.22),
		make_pair(-1.52, -0.88),
		make_pair(-1.51, -1.23),
		make_pair(-1.08, -0.73),
		make_pair(-0.22, -0.94),
		make_pair(-0.50, -0.52),
		make_pair(-1.07, -0.73),
		make_pair(-0.36, -1.06),
		make_pair(-0.81, -0.73),
		make_pair(-6.04, -0.29),
		make_pair(-6.08, -0.75),
		make_pair(-6.05, -0.46),
		make_pair(-4.57, -0.78),
		make_pair(-4.70, -0.99),
		make_pair(-4.72, -1.20),
		make_pair(0.06, -0.93),
		make_pair(0.24, -1.10),
		make_pair(0.08, -0.91),
		make_pair(-0.47, -0.89),
		make_pair(-0.49, -0.88),
		make_pair(-0.40, -0.65),
		make_pair(-0.28, -0.58),
		make_pair(-0.28, -0.60),
		make_pair(-0.26, -0.56),
		make_pair(-0.45, -0.46),
		make_pair(-0.72, -0.20),
		make_pair(-0.62, -0.39),
		make_pair(-5.84, -0.25),
		make_pair(-5.12, -0.74),
		make_pair(-4.98, -1.03),
		make_pair(-5.66, -0.45),
		make_pair(-5.52, -0.25),
		make_pair(-5.38, -0.45),
		make_pair(-3.45, -1.21),
		make_pair(-3.75, -0.98),
		make_pair(-3.49, -1.29),
		make_pair(-1.99, -2.67),
		make_pair(-2.09, -2.60),
		make_pair(-2.08, -2.70),
		make_pair(-5.44, -1.54),
		make_pair(-5.44, -1.55),
		make_pair(-5.68, -1.31),
		make_pair(-4.47, -1.21),
		make_pair(-4.54, -1.12),
		make_pair(-4.82, -1.01),
		make_pair(-4.07, -0.68),
		make_pair(-3.98, -1.17),
		make_pair(-4.07, -0.60),
		make_pair(-3.52, -1.19),
		make_pair(-3.53, -1.22),
		make_pair(-3.57, -1.19),
		make_pair(-5.28, -1.71),
		make_pair(-5.13, -1.86),
		make_pair(-6.13, -0.86),
		make_pair(-6.21, -0.78),
		make_pair(-5.24, -1.75),
		make_pair(-6.77, -0.22),
		make_pair(-2.60, -1.48),
		make_pair(-2.51, -1.59),
		make_pair(-2.43, -1.50),
		make_pair(-2.15, -1.40),
		make_pair(-2.11, -1.44),
		make_pair(-1.91, -1.11),
		make_pair(0.11, -0.62),
		make_pair(-0.10, -0.54),
		make_pair(0.08, -0.70),
		make_pair(-0.01, -0.88),
		make_pair(-0.20, -0.72),
		make_pair(-0.05, -1.00),
		make_pair(-0.13, -0.51),
		make_pair(-0.15, -0.37),
		make_pair(-0.19, -0.38),
		make_pair(-0.24, -0.05),
		make_pair(-0.11, -0.24),
		make_pair(-0.09, -0.25),
		make_pair(-2.57, -1.76),
		make_pair(-1.73, -2.14),
		make_pair(-1.79, -2.28),
		make_pair(-1.16, -1.86),
		make_pair(-1.34, -1.12),
		make_pair(-1.28, -0.90),
		make_pair(-4.59, -1.87),
		make_pair(-4.53, -1.96),
		make_pair(-5.40, -1.06),
		make_pair(-3.26, -1.98),
		make_pair(-2.57, -2.77),
		make_pair(-2.66, -2.77),
		make_pair(-0.24, -0.79),
		make_pair(-0.37, -0.71),
		make_pair(-0.36, -0.74),
		make_pair(-0.03, -1.46),
		make_pair(-0.09, -1.36),
		make_pair(0.01, -1.56),
		make_pair(-1.51, -0.87),
		make_pair(-1.27, -1.52),
		make_pair(-1.36, -1.40),
		make_pair(-1.68, -0.64),
		make_pair(-1.39, -0.92),
		make_pair(-1.63, -0.97),
		make_pair(-5.44, -1.50),
		make_pair(-5.62, -1.12),
		make_pair(-5.58, -0.96),
		make_pair(-6.44, -0.54),
		make_pair(-6.39, -0.60),
		make_pair(-6.32, -0.67),
	};
}
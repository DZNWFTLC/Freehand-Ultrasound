#include <opencv2/opencv.hpp>
#include <opencv2/aruco.hpp>
#include <opencv2/core/matx.hpp>
#include <vector>
#include <opencv2/opencv.hpp>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "iostream"
#include <vector>
#include <stdio.h>
#include <opencv2/imgcodecs.hpp>
#include <string>
#include <math.h>


using namespace cv;
using namespace std;
using namespace aruco;
const double PI = 3.141592653589793238463;


int main(int argc, char** argv)
{
	int wait_time = 10;
	ofstream outFile;

	outFile.open("datanew4.csv", ios::out); // 打开模式可省略

	Mat distortion_coefficients1;
	Mat camera_matrix1;
	Mat distortion_coefficients2;
	Mat camera_matrix2;

	Ptr<DetectorParameters> parameter = DetectorParameters::create();
	parameter->errorCorrectionRate = 0.8f;
	parameter->adaptiveThreshWinSizeStep = 3;
	parameter->adaptiveThreshWinSizeMax = 16;
	parameter->maxMarkerPerimeterRate = 0.5;
	parameter->minMarkerPerimeterRate = 0.01;
	parameter->cornerRefinementMaxIterations = 50;
	parameter->cornerRefinementWinSize = 10;
	parameter->cornerRefinementMinAccuracy = 0.05;

	FileStorage calibration1("cameraL.yaml", FileStorage::READ);
	calibration1["camera_matrix"] >> camera_matrix1;
	calibration1["cistortion_coefficients"] >> distortion_coefficients1;
	calibration1.release();
	VideoCapture video1("..\\Data\\SyncedWithTimestamps\\EndoscopeImageMemory_0_de_test.avi");
	int Frames = 0;
	int Valid = 0;
	while (video1.grab())
	{
		int num_markers = 0;

		Mat image1, image_distortion, image_copy;
		video1.retrieve(image1);

		//undistort(image, image_distortion, camera_matrix, distortion_coefficients);

		image1.copyTo(image_copy);
		vector<int> ids;
		vector<vector<Point2f> > corners;
		//Ptr<aruco::Dictionary> dictionary = aruco::generateCustomDictionary(6, 3);
		Ptr<Dictionary> dictionary = getPredefinedDictionary(DICT_4X4_50);

		//cout << "T1_2  = " << endl << T1_2 << endl << endl;

		detectMarkers(image1, dictionary, corners, ids, parameter);
		//detectMarkers(image2, dictionary, corners2, ids2, parameter);

		// if at least one marker detected
		vector<Vec3d> rvecs, tvecs;
		Vec3d tvec;
		Mat R;
		if (ids.size() > 0) {
			drawDetectedMarkers(image_copy, corners, ids);
			cv::aruco::estimatePoseSingleMarkers(corners, 0.005, camera_matrix1, distortion_coefficients1, rvecs, tvecs);
			for (int i = 0; i < ids.size(); i++) {
				cout << "ID  = " << endl << ids[i] << endl << endl;
				if (ids[i] == 0) {
					tvec = tvecs[i];
					cout << "Tvec  = " << endl << tvec << endl << endl;
					Rodrigues(rvecs[i], R);
					cout << "Rvec  = " << endl << rvecs[i] << endl << endl;
					cout << "corners  = " << endl << corners[i] << endl << endl;

					cout << "R  = " << endl << R << endl << endl;
					Valid++;
					cv::aruco::drawAxis(image_copy, camera_matrix1, distortion_coefficients1, rvecs[i], tvec, 0.01);
					for (int j = 0; j < 3; j++) {
						outFile << tvec[j] << ',';
					}
					for (int j = 0; j < 3; j++) {
						if (j == 2)
							outFile << rvecs[i][j] << endl;
						else
							outFile << rvecs[i][j] << ',';
					}
					/*
					Vec3f Rrow0 = R.row(0);
					Vec3f Rrow1 = R.row(1);
					Vec3f Rrow2 = R.row(2);
					for (int j = 0; j < 3; j++) {
						outFile << Rrow0[j] << ',';
					}
					for (int j = 0; j < 3; j++) {
						outFile << Rrow1[j] << ',';
					}
					for (int j = 0; j < 3; j++) {
						if (j == 2)
							outFile << Rrow2[j] << endl;
						else
							outFile << Rrow2[j] << ',';
					}*/
				}

			}
		}
		Frames++;
		imshow("Detected markers", image_copy);
		if (Frames == 150)
			break;
		char key = (char)waitKey(wait_time);
		if (key == 27)
			break;
	}
	outFile << Valid << ',' << Frames << endl;
	outFile.close();
	video1.release();
}
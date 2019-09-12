/**********************************************************
Name : Aruco marker tracking
Date : 2019/06/10
By   : Chongyun WANG
Final: 
**********************************************************/
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


using namespace cv;
using namespace std;
using namespace aruco;


int main(int argc, char** argv)
{
	ofstream outFile;
	outFile.open("data.csv", ios::out); // 打开模式可省略

	string imageName("..\\Data\\test3.png"); // by default
	if (argc > 1)
	{
		imageName = argv[1];
	}
	Mat image, image_copy;
	image = imread(imageName.c_str(), IMREAD_COLOR); // Read the file
	if (image.empty())                      // Check for invalid input
	{
		cout << "Could not open or find the image" << std::endl;
		return -1;
	}
	Ptr<Dictionary> dictionary = getPredefinedDictionary(DICT_4X4_50);
	vector<int> ids;
	vector<vector<Point2f> > corners, rejectedCandidates;
	image.copyTo(image_copy);
	detectMarkers(image, dictionary, corners, ids);
	// if at least one marker detected
	Vec3d rvecs, tvecs;
	//estimatePoseSingleMarkers(corners, 0.05, cameraMatrix, distCoeffs, rvecs, tvecs);
	if (ids.size() > 0)
		drawDetectedMarkers(image_copy, corners, ids);
	if (ids.size() > 0) {
		for (int i = 0; i < ids.size(); i++)
			outFile << corners[i][0].x << ',' << corners[i][1].x << ',' << corners[i][2].x << ',' << corners[i][3].x << ',' << corners[i][0].y << ',' << corners[i][1].y << ',' << corners[i][2].y << ',' << corners[i][3].y << endl;
	}
	outFile.close();
	imshow("Detected markers", image_copy);

	//waitKey(0); // Wait for a keystroke in the window
	return 0;
	
	/*
	int wait_time = 10;
	VideoCapture in_video;
	in_video.open(0);
	

	while (in_video.grab())
	{
		Mat cameraMatrix, distCoeffs;
		Mat image, image_copy;
		in_video.retrieve(image);
		image.copyTo(image_copy);
		vector<int> ids;
		vector<vector<Point2f> > corners;
		Ptr<Dictionary> dictionary =
			generateCustomDictionary(6, 3);
		//Ptr<Dictionary> dictionary = getPredefinedDictionary(DICT_4X4_50);
		Mat markerImg;
		aruco::drawMarker(dictionary, 1, 200, markerImg, 1);
			imshow("marker", markerImg);

		detectMarkers(image, dictionary, corners, ids);

		// if at least one marker detected
		Vec3d rvecs, tvecs;
		//estimatePoseSingleMarkers(corners, 0.05, cameraMatrix, distCoeffs, rvecs, tvecs);
		if (ids.size() > 0)
			drawDetectedMarkers(image_copy, corners, ids);

		imshow("Detected markers", image_copy);
		char key = (char)waitKey(wait_time);
		if (key == 27)
			break;
	}
	in_video.release();
	*/
	}
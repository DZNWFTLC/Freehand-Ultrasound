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
#include "Eigen/Eigen"

using namespace cv;
using namespace std;
using namespace aruco;
const double PI = 3.141592653589793238463;
using namespace Eigen;

int main(int argc, char** argv)
{
	int wait_time = 10;
	ofstream outFile1;
	ofstream outFile2;

	outFile1.open("datasinglenew3.csv", ios::out); // 打开模式可省略
	outFile2.open("datastereonew3.csv", ios::out); // 打开模式可省略

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

	float marker_size = 0.005;
	Mat distCoeffsL;
	Mat cameraMatrixL;
	Mat distCoeffsR;
	Mat cameraMatrixR;
	Mat projMatL;
	Mat projMatR;

	vector<int> idsL, idsR, ids;
	vector<vector<Point2f>> cornersL, cornersR, corners;
	//Ptr<Dictionary> dictionary = getPredefinedDictionary(DICT_4X4_50);
	Ptr<Dictionary> dictionary = generateCustomDictionary(6, 3);

	FileStorage calibrationL("cameraL.yaml", FileStorage::READ);
	calibrationL["camera_matrix"] >> cameraMatrixL;
	calibrationL["distortion_coefficients"] >> distCoeffsL;
	calibrationL.release();
	cout << "camera_matrixL  = " << endl << cameraMatrixL << endl << endl;
	cout << "distCoeffsL  = " << endl << distCoeffsL << endl << endl;

	FileStorage calibrationR("cameraR.yaml", FileStorage::READ);
	calibrationR["camera_matrix"] >> cameraMatrixR;
	calibrationR["distortion_coefficients"] >> distCoeffsR;
	calibrationR.release();
	cout << "camera_matrixR  = " << endl << cameraMatrixR << endl << endl;
	cout << "distCoeffsR  = " << endl << distCoeffsR << endl << endl;

	FileStorage Stereo_calibration("extrinsics.yml", FileStorage::READ);
	Stereo_calibration["P1"] >> projMatL;
	Stereo_calibration["P2"] >> projMatR;
	Stereo_calibration.release();
	cout << "ProjectionL  = " << endl << projMatL << endl << endl;
	cout << "ProjectionR  = " << endl << projMatR << endl << endl;

	//VideoCapture videoL("..\\Data\\EndoscopeImageMemory_0_de_test.avi");
	//VideoCapture videoR("..\\Data\\EndoscopeImageMemory_1_de_test.avi");
	VideoCapture videoL("..\\Data\\SyncedWithTimestamps\\EndoscopeImageMemory_0_de_test.avi");
	VideoCapture videoR("..\\Data\\SyncedWithTimestamps\\EndoscopeImageMemory_1_de_test.avi");
	int Frames = 0;
	int Valid = 0;

	while (videoL.grab() && videoR.grab())
	{
		int num_markers = 0;

		Mat imageL, imageR, imageL_un, imageR_un, image_copyL, image_copyR;
		videoL.retrieve(imageL);
		videoR.retrieve(imageR);

		//undistort(image, image_distortion, camera_matrix, distortion_coefficients);
		imageL.copyTo(image_copyL);
		imageR.copyTo(image_copyR);


		Mat T_m = Mat::eye(4, 4, CV_64F);
		Mat T1_m = Mat::eye(4, 4, CV_64F);
		Mat Rot_m, Rot_m2;

		undistort(imageL, imageL_un, cameraMatrixL, distCoeffsL);
		undistort(imageR, imageR_un, cameraMatrixR, distCoeffsR);

		detectMarkers(imageL_un, dictionary, cornersL, idsL, parameter);
		detectMarkers(imageR_un, dictionary, cornersR, idsR, parameter);
		detectMarkers(imageL, dictionary, corners, ids, parameter);

		// if at least one marker detected
		vector<Vec3d> rvecs, tvecs;
		//Mat rvec, tvec;
		Vec3d rvec, tvec;

		if (idsL.size() > 0) {
			for (int i = 0; i < idsL.size(); i++) {
				if (idsL[i] == 4) {
					if (idsR.size() > 0) {
						for (int j = 0; j < idsR.size(); j++) {
							if (idsR[j] == 4) {
								estimatePoseSingleMarkers(corners, marker_size, cameraMatrixL, distCoeffsL, rvecs, tvecs);

								drawDetectedMarkers(imageL, corners, ids);
								imshow("Detected markers", imageL);
								for (int k = 0; k < ids.size(); k++) {
									if (ids[k] == 4) {
										rvec = rvecs[k];
										tvec = tvecs[k];
										cout << "tvec  = " << endl << tvec << endl << endl;
										for (int l = 0; l < 3; l++) {
											outFile1 << tvec[l] << ',';
										}
										for (int l = 0; l < 3; l++) {
											if (l == 2)
												outFile1 << rvec[l] << endl;
											else
												outFile1 << rvec[l] << ',';
										}
										cout << "CornersL  = " << endl << cornersL[i] << endl << endl;
										cout << "CornersR  = " << endl << cornersR[j] << endl << endl;
										Valid++;
										outFile2 << cornersL[i][0].x << ',' << cornersL[i][1].x << ',' << cornersL[i][2].x << ',' << cornersL[i][3].x << ',' << cornersL[i][0].y << ',' << cornersL[i][1].y << ',' << cornersL[i][2].y << ',' << cornersL[i][3].y << ',';
										outFile2 << cornersR[j][0].x << ',' << cornersR[j][1].x << ',' << cornersR[j][2].x << ',' << cornersR[j][3].x << ',' << cornersR[j][0].y << ',' << cornersR[j][1].y << ',' << cornersR[j][2].y << ',' << cornersR[j][3].y << endl;
									}
								}

								/*
								Mat triangCoords4D;
								triangulatePoints(projMatL, projMatR, cornersL[i], cornersR[j], triangCoords4D);
								Vec4f triangCoords1 = triangCoords4D.col(0);
								Vec4f triangCoords2 = triangCoords4D.col(1);
								Vec4f triangCoords3 = triangCoords4D.col(2);
								Vec4f triangCoords4 = triangCoords4D.col(3);
								cout << "triangCoords4D  = " << endl << triangCoords4D << endl << endl;
								Vector3f Coords13D, Coords23D, Coords33D, Coords43D;
								for (unsigned int k = 0; k < 3; k++) {
									Coords13D[k] = triangCoords1[k] / triangCoords1[3];
									Coords23D[k] = triangCoords2[k] / triangCoords2[3];
									Coords33D[k] = triangCoords3[k] / triangCoords3[3];
									Coords43D[k] = triangCoords4[k] / triangCoords4[3];
								}
								cout << "Coord1  = " << endl << Coords13D << endl << endl;
								cout << "Coord2  = " << endl << Coords23D << endl << endl;
								cout << "Coord3  = " << endl << Coords33D << endl << endl;
								cout << "Coord4  = " << endl << Coords43D << endl << endl;


								MatrixXf coord_temp(3, 4), coord_origin(3, 4), coord_refined(3, 3), coord, coordCorners;
								coord_temp.col(0) = Coords13D;
								coord_temp.col(1) = Coords23D;
								coord_temp.col(2) = Coords33D;
								coord_temp.col(3) = Coords43D;
								cout << "coord_temp  = " << endl << coord_temp << endl << endl;

								Vector4f trierror;
								for (unsigned int k = 0; k < 4; k++) {
									trierror[k] = (coord_temp.row(2).sum() - coord_temp(2, k)) / 3;
								}

								int valid;
								MatrixXf::Index minRow, minCol;
								coord_temp.minCoeff(&minRow, &minCol);
								if (abs(trierror[minCol]) > marker_size / 2) {
									int element = 0;
									for (unsigned int k = 0; k < 4; k++) {
										if (k != minCol) {
											coord_refined.col(element) = coord_temp.col(k);
											element++;
										}
										else
											continue;
									}
									coord = coord_refined;
									coordCorners = coord_refined.transpose();
								}
								else {
									coord_origin = coord_temp;
									coord = coord_origin;
									coordCorners = coord_origin.transpose();
								}
								cout << "coord  = " << endl << coord << endl << endl;

								cout << "coordCorners  = " << endl << coordCorners << endl << endl;

								Vector3f centroid(coord.row(0).mean(), coord.row(1).mean(), coord.row(2).mean());
								coord.row(0).array() -= centroid(0); coord.row(1).array() -= centroid(1); coord.row(2).array() -= centroid(2);
								auto svd = coord.jacobiSvd(ComputeThinU | ComputeThinV);
								Vector3f plane_normal = svd.matrixU().rightCols<1>();
								cout << "centroid  = " << endl << centroid << endl << endl;
								cout << "plane_normal  = " << endl << plane_normal << endl << endl;
								Vector3f reference(0.0, 0.0, -1.0);

								Vector3f crossProduct = reference.cross(plane_normal);
								float crossProductNorm = crossProduct.norm();
								Eigen::Vector3f vector_X = (crossProduct / crossProductNorm);

								// Step 2: Find angle (theta)
								float dotProduct = crossProduct.dot(plane_normal);
								float norm_A = crossProduct.norm();
								float norm_B = plane_normal.norm();
								float dotProductOfNorms = norm_A * norm_B;
								float dotProductDividedByDotProductOfNorms = (dotProduct / dotProductOfNorms);
								float thetaAngleRad = acos(dotProductDividedByDotProductOfNorms);

								// Step 3: Construct A, the skew-symmetric matrix corresponding to X
								Matrix3f matrix_A = Matrix3f::Identity();

								matrix_A(0, 0) = 0.0;
								matrix_A(0, 1) = -1.0 * (vector_X(2));
								matrix_A(0, 2) = vector_X(1);
								matrix_A(1, 0) = vector_X(2);
								matrix_A(1, 1) = 0.0;
								matrix_A(1, 2) = -1.0 * (vector_X(0));
								matrix_A(2, 0) = -1.0 * (vector_X(1));
								matrix_A(2, 1) = vector_X(0);
								matrix_A(2, 2) = 0.0;

								// Step 4: Plug and chug.
								Matrix3f IdentityMat = Matrix3f::Identity();
								Matrix3f firstTerm = sin(thetaAngleRad) * matrix_A;
								Matrix3f secondTerm = (1.0 - cos(thetaAngleRad)) * matrix_A * matrix_A;
								Matrix3f matrix_R = IdentityMat + firstTerm + secondTerm;

								// This is the rotation matrix. Finished with the Rodrigues' Rotation Formula implementation.
								std::cout << "matrix_R" << std::endl << matrix_R << std::endl;

								// We copy the rotation matrix into the matrix that will be used for the transformation.
								Matrix4f Transform = Matrix4f::Identity();
								Transform(0, 0) = matrix_R(0, 0);
								Transform(0, 1) = matrix_R(0, 1);
								Transform(0, 2) = matrix_R(0, 2);
								Transform(1, 0) = matrix_R(1, 0);
								Transform(1, 1) = matrix_R(1, 1);
								Transform(1, 2) = matrix_R(1, 2);
								Transform(2, 0) = matrix_R(2, 0);
								Transform(2, 1) = matrix_R(2, 1);
								Transform(2, 2) = matrix_R(2, 2);
								Transform(0, 3) = centroid(0);
								Transform(1, 3) = centroid(1);
								Transform(2, 3) = centroid(2);

								cout << "Transform  = " << endl << Transform << endl << endl;

								for (int k = 0; k < 3; k++) {
									if (k == 2)
										outFile2 << Transform(k, 0) << ',' << Transform(k, 1) << ',' << Transform(k, 2) << ',' << Transform(k, 3) << endl;
									else
										outFile2 << Transform(k, 0) << ',' << Transform(k, 1) << ',' << Transform(k, 2) << ',' << Transform(k, 3) << ',';
								}
								*/
							}
						}
					}
					break;
				}
			}
		}
		Frames++;
		if (Frames == 150)
			break;
		char key = (char)waitKey(10);
		if (key == 27)
			break;
	}
	outFile1 << Valid << ',' << Frames << endl;
	outFile1.close();
	outFile2.close();
	videoL.release();
	videoR.release();
}
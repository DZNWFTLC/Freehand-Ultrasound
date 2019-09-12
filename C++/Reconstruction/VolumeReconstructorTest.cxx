/*=IGSIO=header=begin======================================================
  Program: IGSIO
  Copyright (c) Laboratory for Percutaneous Surgery. All rights reserved.
  See License.txt for details.
========================================================IGSIO=header=end*/

#include "igsioConfigure.h"
#include "igsioTrackedFrame.h"
#include "igsioXmlUtils.h"
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"
#include "vtkIGSIOSequenceIO.h"
#include "vtkIGSIOTrackedFrameList.h"
#include "vtkIGSIOTransformRepository.h"
#include "vtkIGSIOVolumeReconstructor.h"
#include "vtkXMLUtilities.h"
#include "vtksys/CommandLineArguments.hxx"
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "ceres/ceres.h"
#include "glog/logging.h"
#include <opencv2/opencv.hpp>
#include <opencv2/aruco.hpp>
#include <opencv2/core/matx.hpp>
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
#include <opencv2/calib3d.hpp>

using namespace cv;
using namespace std;
using namespace aruco;
using namespace Eigen;

const double PI = 3.141592653589793238463;
extern Mat imageL;



int waitTime = 0;
//string dictionaryString;
int iCorrectionRate = 0;

//void cvTackBarEvents(int pos, void*);
int ithresh = 0;

static bool readCameraParameters(string filename, Mat& camMatrix, Mat& distCoeffs) {
	FileStorage fs(filename, FileStorage::READ);
	if (!fs.isOpened())
		return false;
	fs["camera_matrix"] >> camMatrix;
	fs["distortion_coefficients"] >> distCoeffs;
	return true;
}
/*
template<class Vec3f>
pair<Vec3f, Vec3f> best_plane_from_points(const vector<Vec3f>& c)
{
	// copy coordinates to  matrix in Eigen format
	size_t num_atoms = c.size();
	Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic > coord(num_atoms,3);
	for (size_t i = 0; i < num_atoms; ++i) coord.col(i) = c[i];

	// calculate centroid
	Vec3f centroid(coord.row(0).mean(), coord.row(1).mean(), coord.row(2).mean());

	// subtract centroid
	coord.row(0).array() -= centroid(0); coord.row(1).array() -= centroid(1); coord.row(2).array() -= centroid(2);

	// we only need the left-singular matrix here
	//  http://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
	auto svd = coord.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
	Vec3f plane_normal = svd.matrixU().rightCols<1>();
	return std::make_pair(centroid, plane_normal);
}*/

//frames start= 464 end =900


int main(int argc, char** argv) {

	//_________________ Initialize Aruco __________________________//

		// Endoscopic
		// read the input image
	/*
		ofstream outFile;
	outFile.open("data.csv", ios::out); // 打开模式可省略
		if (ids.size() > 0) {
		for (int i = 0; i < ids.size(); i++)
			for (int j = 0; j < 4; j++)
				if (j == 3)
					outFile << corners[i][j].x << ',' << corners[i][j].y << endl;
				else
					outFile << corners[i][j].x << ',' << corners[i][j].y << ',';
	}
	outFile.close();


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

	vector<int> idsL, idsR;
	vector<vector<Point2f>> cornersL, cornersR;
	//Ptr<Dictionary> dictionary = getPredefinedDictionary(DICT_4X4_50);
	Ptr<Dictionary> dictionary = generateCustomDictionary(6, 3);

	FileStorage calibrationL("..\\..\\..\\..\\Source\\VolumeReconstruction\\Testing\\Calibration\\cameraL.yaml", FileStorage::READ);
	calibrationL["camera_matrix"] >> cameraMatrixL;
	calibrationL["distortion_coefficients"] >> distCoeffsL;
	calibrationL.release();
	cout << "camera_matrixL  = " << endl << cameraMatrixL << endl << endl;
	cout << "distCoeffsL  = " << endl << distCoeffsL << endl << endl;

	FileStorage calibrationR("..\\..\\..\\..\\Source\\VolumeReconstruction\\Testing\\Calibration\\cameraR.yaml", FileStorage::READ);
	calibrationR["camera_matrix"] >> cameraMatrixR;
	calibrationR["distortion_coefficients"] >> distCoeffsR;
	calibrationR.release();
	cout << "camera_matrixR  = " << endl << cameraMatrixR << endl << endl;
	cout << "distCoeffsR  = " << endl << distCoeffsR << endl << endl;

	FileStorage Stereo_calibration("..\\..\\..\\..\\Source\\VolumeReconstruction\\Testing\\Calibration\\extrinsics.yml", FileStorage::READ);
	Stereo_calibration["P1"] >> projMatL;
	Stereo_calibration["P2"] >> projMatR;
	Stereo_calibration.release();
	cout << "ProjectionL  = " << endl << projMatL << endl << endl;
	cout << "ProjectionR  = " << endl << projMatR << endl << endl;

	VideoCapture videoL("..\\..\\..\\..\\Source\\VolumeReconstruction\\Testing\\Data\\SyncedWithTimestamps\\EndoscopeImageMemory_0_de_test.avi");
	VideoCapture videoR("..\\..\\..\\..\\Source\\VolumeReconstruction\\Testing\\Data\\SyncedWithTimestamps\\EndoscopeImageMemory_1_de_test.avi");

	int frameNum = 0;
	int validNum = 0;

	while (videoL.grab() && videoR.grab())
	{
		frameNum++;

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
								estimatePoseSingleMarkers(cornersL, marker_size, cameraMatrixL, distCoeffsL, rvecs, tvecs);
								
								drawDetectedMarkers(imageL, cornersL, idsL);
								imshow("Detected markers", imageL);

								rvec = rvecs[i];
								tvec = tvecs[i];

								cout << "CornersL  = " << endl << cornersL[i] << endl << endl;
								cout << "CornersR  = " << endl << cornersR[j] << endl << endl;


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

								validNum++;
								cout << "Valid Frames  = " << endl << validNum << endl << endl;
								cout << "Frames  = " << endl << frameNum << endl << endl;

								MatrixXf coord_temp(3,4), coord_origin(3, 4), coord_refined(3, 3), coord, coordCorners;
								coord_temp.col(0) = Coords13D;
								coord_temp.col(1) = Coords23D;
								coord_temp.col(2) = Coords33D;
								coord_temp.col(3) = Coords43D;
								cout << "coord_temp  = " << endl << coord_temp << endl << endl;

								Vector4f trierror;
								for (unsigned int k = 0; k < 4; k++) {
									trierror[k] = (coord_temp.row(2).sum() - coord_temp(2, k))/3;
								}

								int valid;
								MatrixXf::Index minRow, minCol;
								coord_temp.minCoeff(&minRow, &minCol);
								if (abs(trierror[minCol]) > marker_size / 2) {
									valid = false;
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
									valid = true;
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

								MatrixXf ReferenceMarker(4, 4), TransformedMarker(4, 4);
								Vector4f ReferenceCorner;

								ReferenceMarker << 0.0025, 0.0025, 0, 1, -0.0025, 0.0025, 0, 1, -0.0025, -0.0025, 0, 1, 0.0025, -0.0025, 0, 1;
								for (unsigned int k = 0; k < 4; k++) {
									ReferenceCorner = ReferenceMarker.row(k);
									TransformedMarker.row(k) = Transform *ReferenceCorner;
								}
								MatrixXf PointRefMarker(4, 3);
								cout << "TransformedMarker  = " << endl << TransformedMarker << endl << endl;

								PointRefMarker = TransformedMarker.block(0, 0, 4, 3); 
								cout << "PointRefMarker  = " << endl << PointRefMarker << endl << endl;

								cout << "coordCorners  = " << endl << coordCorners << endl << endl;
								MatrixXf T;

								//ICP_OUT icp_result;

								//icp_result = icp(PointRefMarker, coordCorners, 50, 0.0001);
								//int iter = icp_result.iter;
								//T = icp_result.trans;
								
								if (valid) {
									T = best_fit_transform(PointRefMarker, coordCorners);
									cout << "ICP Transform  = " << endl << T << endl << endl;
								}
								
								//Mat imagePoints;
								//projectPoints(points, rvec, tvec, cameraMatrixL, distCoeffsL, imagePoints);
								//cout << "ImagePoints  = " << endl << imagePoints << endl << endl;
								//imshow("Detect frames", image_copyL);
								
							}
						}
					}
					break;
				}
			}
		}
		char key = (char)waitKey(10);
		if (key == 27)
			break;
	}
	videoL.release();
	videoR.release();
	
	*/


	



	//ofstream outFile;
	//outFile.open("data.csv", ios::out); // 打开模式可省略

	cv::Mat InImage, InImageCopy;
	// Open input and read image
	VideoCapture vreader("..\\..\\..\\..\\Source\\VolumeReconstruction\\Testing\\Data\\datathesis\\EndoscopeImageMemory_0_de_test.avi");
	if (vreader.isOpened())
		vreader >> InImage;
	else
	{
		cerr << "Could not open input" << endl;
		return -1;
	}

	// Ultrasound
	// read the input image
	cv::Mat InImage_Ultra, hist_equalized_Ultraimage;
	// Open input and read image
	VideoCapture vreader_Ultra("..\\..\\..\\..\\Source\\VolumeReconstruction\\Testing\\Data\\datathesis\\RenderedImageMemory_0_test.avi");
	if (vreader_Ultra.isOpened())
		vreader_Ultra >> InImage_Ultra;
	else
	{
		cerr << "Could not open Ultra input" << endl;
		return -1;
	}



	/////////////// End of aruco initialization ///////////////////
	char key = 0;
	//Initialize OpenCV//

	Ptr<aruco::DetectorParameters> detectorParams = aruco::DetectorParameters::create();
	detectorParams->errorCorrectionRate = 0.95f;
	bool showRejected = false;
	bool estimatePose = true;
	float markerLength = 0.005;

	//Ptr<aruco::Dictionary> dictionary = aruco::getPredefinedDictionary(aruco::PREDEFINED_DICTIONARY_NAME::DICT_4X4_50);
	Ptr<Dictionary> dictionary = generateCustomDictionary(6, 3);
	Mat camMatrix, distCoeffs;
	if (true) {
		bool readOk = readCameraParameters("..\\..\\..\\..\\Source\\VolumeReconstruction\\Testing\\Calibration\\cameraL.yaml", camMatrix, distCoeffs);
		if (!readOk) {
			cerr << "Invalid camera 1 file" << endl;
			return 0;
		}
	}


	// End OpenCV initialization //


	//_____________Initialization of IGSIO_______________________//
	bool printHelp(false);
	// Parse command line arguments.

	std::string inputImgSeqFileName;// "C:/Users/Charalampos/Desktop/Elbow_Volume/ElbowUltrasoundSweep.mha";
	std::string inputConfigFileName = "..\\..\\..\\..\\Source\\VolumeReconstruction\\Testing\\Config\\Kidney-VolRec.xml";
	std::string outputVolumeFileName = "Volume_test2.mha";
	std::string outputVolumeAlphaFileNameDeprecated;
	std::string outputVolumeAccumulationFileName;
	std::string outputFrameFileName;
	std::string importanceMaskFileName;
	std::string inputImageToReferenceTransformName = "ImageToReference";


	int verboseLevel = vtkIGSIOLogger::LOG_LEVEL_UNDEFINED;

	bool disableCompression = false;
	vtkIGSIOLogger::Instance()->SetLogLevel(verboseLevel);

	vtkSmartPointer<vtkIGSIOVolumeReconstructor> reconstructor = vtkSmartPointer<vtkIGSIOVolumeReconstructor>::New();

	if (!importanceMaskFileName.empty())
	{
		reconstructor->SetImportanceMaskFilename(importanceMaskFileName);
	}


	LOG_INFO("Reading configuration file:" << inputConfigFileName);
	vtkSmartPointer<vtkXMLDataElement> configRootElement = vtkSmartPointer<vtkXMLDataElement>::New();
	if (igsioXmlUtils::ReadDeviceSetConfigurationFromFile(configRootElement, inputConfigFileName.c_str()) == IGSIO_FAIL)
	{
		LOG_ERROR("Unable to read configuration from file " << inputConfigFileName.c_str());
		//return EXIT_FAILURE;
	}

	if (reconstructor->ReadConfiguration(configRootElement) != IGSIO_SUCCESS)
	{
		LOG_ERROR("Failed to read configuration from " << inputConfigFileName.c_str());
		//return EXIT_FAILURE;
	}

	vtkSmartPointer<vtkIGSIOTransformRepository> transformRepository = vtkSmartPointer<vtkIGSIOTransformRepository>::New();
	if (configRootElement->FindNestedElementWithName("CoordinateDefinitions") != NULL)
	{
		if (transformRepository->ReadConfiguration(configRootElement) != IGSIO_SUCCESS)
		{
			LOG_ERROR("Failed to read transforms from CoordinateDefinitions");
			//return EXIT_FAILURE;
		}
	}
	else
	{
		LOG_DEBUG("No transforms were found in CoordinateDefinitions. Only the transforms defined in the input image will be available.");
	}

	// Print calibration transform
	std::ostringstream osTransformRepo;
	transformRepository->Print(osTransformRepo);
	LOG_DEBUG("Transform repository: \n" << osTransformRepo.str());

	// Read image sequence
	LOG_INFO("Reading image sequence " << inputImgSeqFileName);
	vtkSmartPointer<vtkIGSIOTrackedFrameList> trackedFrameList = vtkSmartPointer<vtkIGSIOTrackedFrameList>::New();

	double timestamp = 0;
	int images_count = 0;
	int used_frames = 200;

	//Transformations:
	float aray[80][16]; // here data will be saved, data from txt file or streaming-tracking
	int i_transformations = 0;
	///////////// End of IGSIO initialization ///////////////////

	bool detection_happened = false;
	int frames_counter = 0;

	while (vreader.grab() && vreader_Ultra.grab()&& frames_counter < used_frames)
	{
		/////////////// Start Aruco tracking and obtaining the matrices ////////////////////////////////

		vreader_Ultra.retrieve(InImage_Ultra);
		//cvtColor(InImage_Ultra, InImage_Ultra, COLOR_RGB2GRAY);
		//equalizeHist(InImage_Ultra, hist_equalized_Ultraimage);
		vreader.retrieve(InImage);

		InImage.copyTo(InImageCopy);

		// Let's detect with OpenCV!!!

		vector< int > ids;
		vector< vector< Point2f > > corners, rejected;
		vector< Vec3d > rvecs, tvecs;
		//Mat tvecs;
		// detect markers and estimate pose
		aruco::detectMarkers(InImage, dictionary, corners, ids, detectorParams, rejected);
		Mat R;
		Vec3f tvec;


			if (ids.size() > 0)
			{
				for (int i = 0; i < ids.size(); i++) {
					cout << "ID  = " << endl << ids[i] << endl << endl;
					if (ids[i] == 4) {
						detection_happened = true;
						aruco::estimatePoseSingleMarkers(corners, markerLength, camMatrix, distCoeffs, rvecs,
							tvecs);
						Rodrigues(rvecs[i], R);
						tvec = tvecs[i];
						cout << "tvecs  = " << endl << tvec << endl << endl;
						/*
						for (int j = 0; j < 3; j++){
								outFile << tvec[j] << ',' ;
					}
						for (int j = 0; j < 3; j++){
							if (j == 2)
							outFile << rvecs[i][j] << endl;
							else
							outFile << rvecs[i][j] << ',';
						}*/
						Vec3d t4(-0.025, 0.025, -0.056);
						Vec3d t0(-0.275, 0.025, -0.056);
					}

				}
			}


			// draw results



			//////////// End of Aruco Tracking /////////////////////////

			/////////// Obtain the coressponding Ultrasound frame /////
			///Images
			igsioVideoFrame& video = igsioVideoFrame();

			if (detection_happened) {
				char name[100];
				if (images_count < 10) {
					snprintf(name, sizeof name, "..\\..\\..\\..\\Source\\VolumeReconstruction\\Testing\\ModelData\\Kidney000%d.bmp", images_count);
				}
				else if (images_count < 100) {
					snprintf(name, sizeof name, "..\\..\\..\\..\\Source\\VolumeReconstruction\\Testing\\ModelData\\Kidney00%d.bmp", images_count);


				}
				else if (images_count < 1000) {
					snprintf(name, sizeof name, "..\\..\\..\\..\\Source\\VolumeReconstruction\\Testing\\ModelData\\Kidney0%d.bmp", images_count);
				}
				else {
					snprintf(name, sizeof name, "..\\..\\..\\..\\Source\\VolumeReconstruction\\Testing\\ModelData\\Kidney%d.bmp", images_count);
				}
				images_count++;

				imwrite(name, InImage_Ultra);


				//cv::imshow("in2", InImage_Ultra);
				///////////////////////////////////////////////////////////

				/////////////// Start making frames for IGSIO /////////////

				igsioTrackedFrame trackedframe = igsioTrackedFrame();


				std::string tring = std::string(name);
				const char* fileName = tring.c_str();
				igsioStatus stat = igsioVideoFrame::ReadImageFromFile(video, fileName);
				trackedframe.SetImageData(video);

				///Matrices
				vtkMatrix4x4* matrix4_Probe_Ref = vtkMatrix4x4::New();

				//matrix4_Probe_Ref->SetElement(0, 0, *R.ptr<double>(0, 0));

				for (int k = 0; k < 3; k++) {
					for (int j = 0; j < 3; j++) {
						matrix4_Probe_Ref->SetElement(k, j, *R.ptr<double>(k, j));
					}
				}
				matrix4_Probe_Ref->SetElement(0, 3, tvec[0]);
				matrix4_Probe_Ref->SetElement(1, 3, tvec[1]);
				matrix4_Probe_Ref->SetElement(2, 3, tvec[2]);
				matrix4_Probe_Ref->SetElement(3, 3, 1.);

				matrix4_Probe_Ref->SetElement(3, 0, 0.f);
				matrix4_Probe_Ref->SetElement(3, 1, 0.f);
				matrix4_Probe_Ref->SetElement(3, 2, 0.f);




				const igsioTransformName  nameProbRef = igsioTransformName(std::string("Probe"), std::string("Reference"));
				trackedframe.SetFrameTransform(nameProbRef, matrix4_Probe_Ref);
				trackedframe.SetTimestamp(timestamp);
				trackedframe.SetFrameTransformStatus(nameProbRef, ToolStatus::TOOL_OK);

				trackedFrameList->AddTrackedFrame(&trackedframe);


				timestamp = timestamp + 0.08;
				frames_counter++;

				detection_happened = false;
			}
			/////////// End of frames ////////////////////////////////


		}  
		//outFile.close();

		//////////// Start of reconstruction algorithm /////////////////////

		igsioTransformName imageToReferenceTransformName;
		if (!inputImageToReferenceTransformName.empty())
		{
			// image to reference transform is specified at the command-line
			if (imageToReferenceTransformName.SetTransformName(inputImageToReferenceTransformName.c_str()) != IGSIO_SUCCESS)
			{
				LOG_ERROR("Invalid image to reference transform name: " << inputImageToReferenceTransformName);
				//return EXIT_FAILURE;
			}
			reconstructor->SetImageCoordinateFrame(imageToReferenceTransformName.From());
			reconstructor->SetReferenceCoordinateFrame(imageToReferenceTransformName.To());
		}

		LOG_INFO("Set volume output extent...");
		std::string errorDetail;
		if (reconstructor->SetOutputExtentFromFrameList(trackedFrameList, transformRepository, errorDetail) != IGSIO_SUCCESS)
		{
			LOG_ERROR("Failed to set output extent of volume!");
			//return EXIT_FAILURE;
		}

		LOG_INFO("Reconstruct volume...");
		const int numberOfFrames = trackedFrameList->GetNumberOfTrackedFrames();
		int numberOfFramesAddedToVolume = 0;

		for (int frameIndex = 0; frameIndex < numberOfFrames; frameIndex += reconstructor->GetSkipInterval())
		{
			LOG_DEBUG("Frame: " << frameIndex);
			vtkIGSIOLogger::PrintProgressbar((100.0 * frameIndex) / numberOfFrames);

			igsioTrackedFrame* frame = trackedFrameList->GetTrackedFrame(frameIndex);

			if (transformRepository->SetTransforms(*frame) != IGSIO_SUCCESS)
			{
				LOG_ERROR("Failed to update transform repository with frame #" << frameIndex);
				continue;
			}

			// Insert slice for reconstruction
			bool insertedIntoVolume = false;
			if (reconstructor->AddTrackedFrame(frame, transformRepository, &insertedIntoVolume) != IGSIO_SUCCESS)
			{
				LOG_ERROR("Failed to add tracked frame to volume with frame #" << frameIndex);
				continue;
			}

			if (insertedIntoVolume)
			{
				numberOfFramesAddedToVolume++;
			}

			// Write an ITK image with the image pose in the reference coordinate system
			if (!outputFrameFileName.empty())
			{
				vtkSmartPointer<vtkMatrix4x4> imageToReferenceTransformMatrix = vtkSmartPointer<vtkMatrix4x4>::New();
				if (transformRepository->GetTransform(imageToReferenceTransformName, imageToReferenceTransformMatrix) != IGSIO_SUCCESS)
				{
					std::string strImageToReferenceTransformName;
					imageToReferenceTransformName.GetTransformName(strImageToReferenceTransformName);
					LOG_ERROR("Failed to get transform '" << strImageToReferenceTransformName << "' from transform repository!");
					continue;
				}

				// Print the image to reference transform
				std::ostringstream os;
				imageToReferenceTransformMatrix->Print(os);
				LOG_TRACE("Image to reference transform: \n" << os.str());

				// Insert frame index before the file extension (image.mha => image001.mha)
				std::ostringstream ss;
				size_t found;
				found = outputFrameFileName.find_last_of(".");
				ss << outputFrameFileName.substr(0, found);
				ss.width(3);
				ss.fill('0');
				ss << frameIndex;
				ss << outputFrameFileName.substr(found);

				//igsioCommon::WriteToFile(frame, ss.str(), imageToReferenceTransformMatrix);
			}
		}

		vtkIGSIOLogger::PrintProgressbar(100);

		trackedFrameList->Clear();

		LOG_INFO("Number of frames added to the volume: " << numberOfFramesAddedToVolume << " out of " << numberOfFrames);

		LOG_INFO("Saving volume to file...");
		reconstructor->SaveReconstructedVolumeToFile(outputVolumeFileName, false, !disableCompression);


		///////////// END OF RECONSTRUCTION ////////////////////////////////

		///////////// START OF RENDERING //////////////////////////////////

		//////////// END OF RENDERING	//////////////////////////////////
		return EXIT_SUCCESS;
}

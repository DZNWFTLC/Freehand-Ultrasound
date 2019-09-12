
// Include for IGSIO //
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
#include "igsioVideoFrame.h"
#include "vtkDataArray.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include <vtkCamera.h>
#include <vtkColorTransferFunction.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkMetaImageReader.h>
#include <vtkNamedColors.h>
#include <vtkPiecewiseFunction.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkVolume.h>
#include <vtkVolumeProperty.h>



// Include for ARUCO //
#include <string>
#include <windows.h>

//#include "aruco.h"
//#include "cvdrawingutils.h"
#include <iostream>
//#include <opencv2/highgui/highgui.hpp>

//Include for OpenCV
#include <opencv2/aruco.hpp>
#include <opencv2/highgui.hpp>
//#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/calib3d.hpp>


#define INF 10000 


using namespace std;
using namespace cv;
//using namespace aruco;

int waitTime = 0;
//string dictionaryString;
int iCorrectionRate = 0;

//void cvTackBarEvents(int pos, void*);
int ithresh = 0;




string ExePath() {
	char buffer[MAX_PATH];
	GetModuleFileName(NULL, buffer, MAX_PATH);
	string::size_type pos = string(buffer).find_last_of("\\/");
	return string(buffer).substr(0, pos);
}


static bool readCameraParameters(string filename, Mat &camMatrix, Mat &distCoeffs) {
	FileStorage fs(filename, FileStorage::READ);
	if (!fs.isOpened())
		return false;
	fs["camera_matrix"] >> camMatrix;
	fs["distortion_coefficients"] >> distCoeffs;
	return true;
}


char* window_name = "Threshold Demo";

char* trackbar_type = "Type: \n 0: Binary \n 1: Binary Inverted \n 2: Truncate \n 3: To Zero \n 4: To Zero Inverted";
char* trackbar_value = "Value";

int threshold_value = 0;
int threshold_type = 3;;
int const max_value = 255;
int const max_type = 4;
int const max_BINARY_value = 255;


Mat new_cut_im , dst;

const int w = 500;
int levels = 3;

vector<vector<Point> > contours;
vector<Vec4i> hierarchy;



int thresh = 100;
RNG rng(12345);
Mat src_gray;
Mat canny_output;
Mat drawing;


Mat drawing2 = Mat::zeros(230,170, CV_8UC3);

// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr' 
bool onSegment(Point p, Point q, Point r)
{
	if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) &&
		q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))
		return true;
	return false;
}

// To find orientation of ordered triplet (p, q, r). 
// The function returns following values 
// 0 --> p, q and r are colinear 
// 1 --> Clockwise 
// 2 --> Counterclockwise 
int orientation(Point p, Point q, Point r)
{
	int val = (q.y - p.y) * (r.x - q.x) -
		(q.x - p.x) * (r.y - q.y);

	if (val == 0) return 0;  // colinear 
	return (val > 0) ? 1 : 2; // clock or counterclock wise 
}

// The function that returns true if line segment 'p1q1' 
// and 'p2q2' intersect. 
bool doIntersect(Point p1, Point q1, Point p2, Point q2)
{
	// Find the four orientations needed for general and 
	// special cases 
	int o1 = orientation(p1, q1, p2);
	int o2 = orientation(p1, q1, q2);
	int o3 = orientation(p2, q2, p1);
	int o4 = orientation(p2, q2, q1);

	// General case 
	if (o1 != o2 && o3 != o4)
		return true;

	// Special Cases 
	// p1, q1 and p2 are colinear and p2 lies on segment p1q1 
	if (o1 == 0 && onSegment(p1, p2, q1)) return true;

	// p1, q1 and p2 are colinear and q2 lies on segment p1q1 
	if (o2 == 0 && onSegment(p1, q2, q1)) return true;

	// p2, q2 and p1 are colinear and p1 lies on segment p2q2 
	if (o3 == 0 && onSegment(p2, p1, q2)) return true;

	// p2, q2 and q1 are colinear and q1 lies on segment p2q2 
	if (o4 == 0 && onSegment(p2, q1, q2)) return true;

	return false; // Doesn't fall in any of the above cases 
}

// Returns true if the point p lies inside the polygon[] with n vertices 
bool isInside(vector<Point> polygon, int n, Point p)
{
	// There must be at least 3 vertices in polygon[] 
	if (n < 3)  return false;
	
	// Create a point for line segment from p to infinite 
	Point extreme = { INF, p.y-50 };
	circle(drawing2, Point(extreme.x, extreme.y), 2, Scalar(0, 0, 255), FILLED, LINE_AA);

	// Count intersections of the above line with sides of polygon 
	int count = 0, i = 0;
	do
	{
		//int next = (i + 1) ;

		// Check if the line segment from 'p' to 'extreme' intersects 
		// with the line segment from 'polygon[i]' to 'polygon[next]' 

		if (doIntersect(polygon[i], polygon[i+1], p, extreme))
		{
			// If the point 'p' is colinear with line segment 'i-next', 
			// then check if it lies on segment. If it lies, return true, 
			// otherwise false 
			if (orientation(polygon[i], p, polygon[i+1]) == 0)
				return onSegment(polygon[i], p, polygon[i+1]);

			count++;
		}
		i = (i + 1) % n;
	} while (i !=0);

	// Return true if count is odd, false otherwise 
	return count & 1;  // Same as (count%2 == 1) 
}


int main(int argc, char** argv)
{
	//_________________ Initialize Aruco __________________________//
	
	cout << "my directory is " << ExePath() << "\n";

		//aruco::CameraParameters CamParam;

		// Endoscopic
		// read the input image
		cv::Mat InImage ,InImageCopy;
		// Open input and read image
		// Phantom 1
		//VideoCapture vreader("C:/Users/Charalampos/Desktop/Ucl/Data_13_6/Log_D2P282649_2019.06.12_16.38.12/part0003/SyncedWithTimestamps/EndoscopeImageMemory_0_crop_3.avi");
		// Phantom 2 Cropped: 0-2436
		VideoCapture vreader("C:/IGSIO-master/IGSIO-master/data/EndoscopeImageMemory_Crop_100_frames.avi");

		if (vreader.isOpened())
			vreader >> InImage;
		else
		{
			cerr << "Could not open input" << endl;
			return -1;
		}
		
		// Ultrasound
		// read the input image
		cv::Mat InImage_Ultra;
		// Open input and read image
		// Phantom 1
		//VideoCapture vreader_Ultra("C:/Users/Charalampos/Desktop/Ucl/Data_13_6/Log_D2P282649_2019.06.12_16.38.12/part0003/SyncedWithTimestamps/RenderedImageMemory_0_crop_3.avi");
		// Phantom 2
		VideoCapture vreader_Ultra("C:/IGSIO-master/IGSIO-master/data/RenderedImageMemory_0_Crop_100_frames.avi");

		if (vreader_Ultra.isOpened())
			vreader_Ultra >> InImage_Ultra;
		else
		{
			cerr << "Could not open Ultra input" << endl;
			return -1;
		}

		// read camera parameters if specifed
		/*CamParam.readFromXMLFile("C:/opencv-master/opencv-master/samples/cpp/tutorial_code/calib3d/camera_calibration/out_camera_data_15_6.xml");
		CamParam.resize(InImage.size());
		// read marker size if specified (default value -1)
		float MarkerSize = 0.005;
		// Create the detector
		MarkerDetector MDetector;
		dictionaryString = cml("-d", "C:/aruco-master/aruco-master/aruco-1.3.0-testsdata/hrm_dictionaries/d4x4_100.dict");
		MDetector.setDictionary(dictionaryString
			, 0.8f);
		MDetector.setThresholdParamRange(2, 0);

		MDetector.setThresholdMethod(MarkerDetector::ThresholdMethods::FIXED_THRES);
		MDetector.setThresholdParams(70, 7);

		std::map<uint32_t, MarkerPoseTracker>
			MTracker;  // use a map so that for each id, we use a different pose tracker
					   // Set the dictionary you want to work with, if you included option -d in command line
					   // se
		char key = 0;*/

		/////////////// End of aruco initialization ///////////////////
		char key = 0; 
		//Initialize OpenCV//

		Ptr<aruco::DetectorParameters> detectorParams = aruco::DetectorParameters::create();
		detectorParams->errorCorrectionRate = 0.8f;
		bool showRejected = false;
		bool estimatePose = true;
		float markerLength = 0.005;

		//Ptr<aruco::Dictionary> dictionary =
		//	aruco::getPredefinedDictionary(aruco::PREDEFINED_DICTIONARY_NAME::DICT_ARUCO_ORIGINAL);

		Ptr<aruco::Dictionary> dictionary =
			aruco::generateCustomDictionary(6, 3);

		Mat camMatrix, distCoeffs;
		if (true) {
			bool readOk = readCameraParameters("C:/opencv-master/opencv-master/samples/cpp/tutorial_code/calib3d/camera_calibration/Camera_1_calibration.xml", camMatrix, distCoeffs);
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
		std::string inputConfigFileName = "C:/Users/Charalampos/Desktop/Elbow_Volume/Kidney-VolRec.xml";
		std::string outputVolumeFileName = "C:/IGSIO-master/IGSIO-master/data/Volume_xaris_21_8.mha";
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
		
		double timestamp = 232.542071;
		int images_count = 0;

		//Transformations:
		float aray[80][16]; // here data will be saved, data from txt file or streaming-tracking
		int i_transformations = 0;
		///////////// End of IGSIO initialization ///////////////////

		bool detection_happened = false;
		int frames_counter = 0;
		int used_frames = 30;
		int first_key = 0;
		std::vector<KeyPoint> first_keypoints;
		int point_count = 0;

		bool detectTumor = false;

		do
		{

			/////////////// Start Aruco tracking and obtaining the matrices ////////////////////////////////
			
			vreader_Ultra.retrieve(InImage_Ultra);

			vreader.retrieve(InImage);

			// Let's detect with ARUCO!!!
			/*vector<Marker> Markers = MDetector.detect(InImage);
			for (auto& marker : Markers)  // for each marker
			{	
				if (marker.id == 4) {
					detection_happened=MTracker[marker.id].estimatePose(marker, CamParam,
						MarkerSize);  // call its tracker and estimate the pose
				}
				else {
					MTracker[marker.id].estimatePose(marker, CamParam,
						MarkerSize);  // call its tracker and estimate the pose
				}
			}
		
			for (unsigned int i = 0; i < Markers.size(); i++)
			{
				//cout << Markers[i] << endl;	// print translation and rotation matrices
				Markers[i].draw(InImage, Scalar(0, 0, 255), 2);
			}
			// draw a 3d cube in each marker if there is 3d info
			if (CamParam.isValid() && MarkerSize != -1)
			{
				for (unsigned int i = 0; i < Markers.size(); i++)
				{
					CvDrawingUtils::draw3dCube(InImage, Markers[i], CamParam);
					CvDrawingUtils::draw3dAxis(InImage, Markers[i], CamParam);
				}
			}
			// show input with augmented information
			cv::namedWindow("in", 1);
			cv::imshow("in", InImage);
			//cv::imshow("in", InImage_Ultra);

			
			*/
			key = cv::waitKey(waitTime);  // wait for key to be pressed
			if (key == 'r')
				waitTime = waitTime == 0 ? 1 : 0;
			// Let's detect with OpenCV!!!

			vector< int > ids;
			vector< vector< Point2f > > corners, rejected;
			vector< Vec3d > rvecs, tvecs;
			//Mat tvecs;
			// detect markers and estimate pose
			aruco::detectMarkers(InImage, dictionary, corners, ids, detectorParams, rejected);
			
			if (estimatePose && ids.size() > 0 ) {

				if (ids[0] == 4)
				{
					detection_happened = true;
					aruco::estimatePoseSingleMarkers(corners, markerLength, camMatrix, distCoeffs, rvecs,
						tvecs);
					aruco::drawAxis(InImage, camMatrix, distCoeffs, rvecs[0], tvecs[0],
						markerLength );
					
				}
				
			}


			// draw results
			InImage.copyTo(InImageCopy);
			if (ids.size() > 0) {
				aruco::drawDetectedMarkers(InImageCopy, corners, ids);

				if (estimatePose && detection_happened) {
					for (unsigned int i = 0; i < ids.size(); i++) {
						if (ids[0] == 4) {
							//aruco::drawAxis(InImageCopy, camMatrix, distCoeffs,rvecs, tvecs,
							//	markerLength * 0.5f);
						}
					}
				}
			}

			imshow("out", InImageCopy);


			//////////// End of Aruco Tracking /////////////////////////

			/////////// Obtain the coressponding Ultrasound frame /////
			///Images
			igsioVideoFrame& video = igsioVideoFrame();
			if (detection_happened) {
				char name[100];
				if (images_count < 10) {
					snprintf(name, sizeof name, "C:/IGSIO-master/IGSIO-master/build/Debug/Kidney_21_8/Kidney000%d.bmp", images_count);
				}
				else if (images_count<100){
					snprintf(name, sizeof name, "C:/IGSIO-master/IGSIO-master/build/Debug/Kidney_21_8/Kidney00%d.bmp", images_count);


				}
				else if (images_count<1000) {
					snprintf(name, sizeof name, "C:/IGSIO-master/IGSIO-master/build/Debug/Kidney_21_8/Kidney0%d.bmp", images_count);


				}
				images_count++;
		
				//imwrite(name, InImage_Ultra);
				

				//cv::imshow("in2", InImage_Ultra);
				///////////////////////////////////////////////////////////

				/////////////// Start making frames for IGSIO /////////////

				igsioTrackedFrame trackedframe = igsioTrackedFrame();

				//MTracker[4].getRTMatrix().empty();
				//Check if tracking works

				
				// Transformation matrix  ARUCO
				//Mat mm = MTracker[4].getRTMatrix();
				
				// Transformation matrix  OpenCV
				Mat mm, R;
				Mat gray;
				Rodrigues(rvecs[0], R);
				Rect RecUltra(614, 72, 1180, 734);

				Mat onlyUltra = InImage_Ultra(RecUltra);

				//---------------------------------------------------------------------------------------------//
				//---------------------------------------Image Segmentation-----------------------------------//
				//--------------------------------------------------------------------------------------------//

				if (detectTumor) {
					
					Rect Rec(1150, 280, 230, 210);
					Rect Rec2(1160, 285, 230, 170);
					//Mat Roi = InImage_Ultra(Rec);

					Mat image1 = InImage_Ultra;// = imread("C:/IGSIO-master/IGSIO-master/build/Debug/Kidney_17_8/Kidney0024.bmp");
					Mat image = image1(Rec2);
					cvtColor(image, gray, cv::COLOR_BGR2GRAY);
					InImage_Ultra = gray;


					// medfilt2
					for (int i = 1; i < 25; i = i + 2)
					{
						medianBlur(gray, dst, i);

					}

					Mat dst1, dst2, dst3, dst4;

					// Not used ( adapthisteq)
					cv::Ptr<cv::CLAHE> clahe = cv::createCLAHE();
					clahe->apply(dst, dst1);

					// imsharpen
					Mat blurred; double sigma = 20, threshold_ = 5, amount = 2;
					GaussianBlur(dst, blurred, Size(), sigma, sigma);
					Mat lowContrastMask = abs(dst - blurred) < threshold_;
					Mat sharpened = dst*(1 + amount) + blurred*(-amount);
					dst.copyTo(sharpened, lowContrastMask);

					// A_sharp - A3
					dst3 = sharpened - dst;

					// imbinarize
					threshold(dst3, dst4, 0, 255, cv::THRESH_BINARY);


					// Set up the detector with default parameters.
					SimpleBlobDetector::Params params;


					params.blobColor = 255;
					// Change thresholds
					//params.minThreshold = 10;
					//params.maxThreshold = 200;

					// Filter by Area.
					params.filterByArea = true;
					params.minArea = 2000;
					params.maxArea = 100000;

					// Filter by Circularity
					params.filterByCircularity = true;
					params.minCircularity = 0.1;

					// Filter by Convexity
					params.filterByConvexity = true;
					params.minConvexity = 0.15;

					// Filter by Inertia
					params.filterByInertia = true;
					params.minInertiaRatio = 0.01;
					// Detect blobs.
					std::vector<KeyPoint> keypoints;
					std::vector<KeyPoint> interest_keypoints;


					cv::Ptr<cv::SimpleBlobDetector> detector = cv::SimpleBlobDetector::create(params);
					detector->detect(dst4, keypoints);

					// Draw detected blobs as red circles.
					// DrawMatchesFlags::DRAW_RICH_KEYPOINTS flag ensures the size of the circle corresponds to the size of blob
					Mat im_with_keypoints;
					if (keypoints.size() > 0 && first_key == 0) {
						first_keypoints = keypoints;
						first_key = 1;
					}

					if (keypoints.size() > 0) {
						for (int i = 0; i < keypoints.size(); i++) {


							if ((keypoints.at(i).pt.x < first_keypoints.at(0).pt.x + 40) && (keypoints.at(i).pt.x > first_keypoints.at(0).pt.x - 40) && (keypoints.at(i).pt.y < first_keypoints.at(0).pt.y + 40) && (keypoints.at(i).pt.y > first_keypoints.at(0).pt.y - 40))
								interest_keypoints.push_back(keypoints.at(i));

						}
					}

					if (interest_keypoints.size() > 0) {
						drawKeypoints(gray, interest_keypoints, im_with_keypoints, Scalar(0, 0, 255), DrawMatchesFlags::DRAW_RICH_KEYPOINTS);


					}
					else {
						SimpleBlobDetector::Params params2;


						params2.blobColor = 0;
						// Change thresholds
						//params.minThreshold = 10;
						//params.maxThreshold = 200;

						// Filter by Area.
						params2.filterByArea = true;
						params2.minArea = 2000;
						params2.maxArea = 100000;

						// Filter by Circularity
						params2.filterByCircularity = true;
						params2.minCircularity = 0.1;

						// Filter by Convexity
						params2.filterByConvexity = true;
						params2.minConvexity = 0.15;

						// Filter by Inertia
						params2.filterByInertia = true;
						params2.minInertiaRatio = 0.01;
						// Detect blobs.


						cv::Ptr<cv::SimpleBlobDetector> detector2 = cv::SimpleBlobDetector::create(params2);
						detector2->detect(dst4, keypoints);

						if (keypoints.size() > 0) {
							for (int i = 0; i < keypoints.size(); i++) {


								if ((keypoints.at(i).pt.x < first_keypoints.at(0).pt.x + 40) && (keypoints.at(i).pt.x > first_keypoints.at(0).pt.x - 40) && (keypoints.at(i).pt.y < first_keypoints.at(0).pt.y + 40) && (keypoints.at(i).pt.y > first_keypoints.at(0).pt.y - 40))
									interest_keypoints.push_back(keypoints.at(i));

							}
						}
						if (interest_keypoints.size() > 0) {
							drawKeypoints(gray, interest_keypoints, im_with_keypoints, Scalar(0, 0, 255), DrawMatchesFlags::DRAW_RICH_KEYPOINTS);


						}


					}
					//imshow("keypoints2", dst4);


					// Show blobs
					//imshow("keypoints", im_with_keypoints);

					vector<vector<Point> > contours2;

					if (interest_keypoints.size() > 0) {

						/// Detect edges using canny
						Canny(dst4, canny_output, 20, 20 * 2, 7);
						/// Find contours
						findContours(canny_output, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, Point(0, 0));
						
						/// Draw contours
						Mat drawing = Mat::zeros(canny_output.size(), CV_8UC3);
						Mat labels;

						std::vector<Point2f> centers;
						point_count = 0;
						Mat points = cv::Mat(contours.size(), 2, CV_32F);

						for (int i = 0; i < contours.size(); i++)
						{
							if (contours[i].at(0).x <interest_keypoints.at(0).pt.x + 90 && contours[i].at(0).x >interest_keypoints.at(0).pt.x - 90 && contours[i].at(0).y <+interest_keypoints.at(0).pt.y + 90 && contours[i].at(0).y >interest_keypoints.at(0).pt.y - 90)
							{
							
								points.at<float>(point_count, 0) = contours[i].at(0).x;
								points.at<float>(point_count, 1) = contours[i].at(0).y;
								
								point_count++;
							}
						}

						/*float sumx = 0;
						float sumy = 0;

						for (int i = 0; i < point_count; i++) {

							sumx = sumx + points.at<float>(i, 0);
							sumy = sumy + points.at<float>(i, 1);

						}
						float meanx = sumx / point_count;
						float meany = sumy / point_count;
						Mat points2 = cv::Mat(point_count, 2, CV_32F);

						int count_2 = 0;
						for (int i = 0; i < point_count; i++) {

							if (points.at<float>(i, 0) < meanx + 30 && points.at<float>(i, 0) > meanx - 30 && points.at<float>(i, 1) < meany + 30 && points.at<float>(i, 1) > meany + 30)
							{
								points2.at<float>(count_2, 0) = points.at<float>(i, 0);
								points2.at<float>(count_2, 1) = points.at<float>(i, 1);

								count_2++;
							}

						}
						for (int i = 0; i < points2.rows; i++)
						{
							int clusterIdx = labels.at<int>(i);
							//Point ipt = points2.at<Point2f>(i);
							Point ipt = Point(points2.at<float>(i, 0), points2.at<float>(i, 1));

							circle(drawing, ipt, 2, Scalar(0, 0, 255), FILLED, LINE_AA);
						}*/



						/*for (int i = 0; i < points2.rows; i++) {
							for (int j = 0; j < 2; j++) {

								points2.at<float>(i, j) = points.at<float>(i, j);

							}

						}*/


						/*	double compactness = kmeans(points2, 2, labels,
								TermCriteria(TermCriteria::EPS + TermCriteria::MAX_ITER, 10, 1.0),
								10, KMEANS_RANDOM_CENTERS, centers);


							for (int i = 0; i < points2.rows; i++)
							{
								int clusterIdx = labels.at<int>(i);
								//Point ipt = points2.at<Point2f>(i);
								Point ipt = Point(points2.at<float>(i,0), points2.at<float>(i, 1));

								circle(drawing, ipt, 2, Scalar(0, 0, 255), FILLED, LINE_AA);
							}
							for (int i = 0; i < (int)centers.size(); ++i)
							{
								Point2f c = centers[i];
								circle(drawing, c, 40, Scalar(0, 255, 180), 1, LINE_AA);
							}
							*/
							//cv::connectedComponents(drawing, labels, 8, CV_32S);
							/// Find the rotated rectangles and ellipses for each contour

						vector<RotatedRect> minRect(contours.size());
						vector<RotatedRect> minEllipse(contours.size());
						RotatedRect tryEllipse;

						for (int i = 0; i < contours.size(); i++)
						{
							minRect[i] = minAreaRect(Mat(contours[i]));
							if (minRect[i].center.x <interest_keypoints.at(0).pt.x + 40 && minRect[i].center.x >interest_keypoints.at(0).pt.x - 40 && minRect[i].center.y <interest_keypoints.at(0).pt.y + 40 && minRect[i].center.y >interest_keypoints.at(0).pt.y - 40) {
								if (contours[i].size() > 5)
								{
									tryEllipse = fitEllipse(Mat(contours[i]));
									if (tryEllipse.center.x <interest_keypoints.at(0).pt.x + 30 && tryEllipse.center.x >interest_keypoints.at(0).pt.x - 30 && tryEllipse.center.y <interest_keypoints.at(0).pt.y + 30 && tryEllipse.center.y >interest_keypoints.at(0).pt.y - 30) {
										minEllipse[i] = tryEllipse;

									}
								}
							}
						}
						vector<vector<Point>>  best_contours(contours.size());
						RotatedRect bestEllipse;
						float max_area = 0;

						for (int i = 0; i < contours.size(); i++) {

							//Size2f wh = minEllipse[i].size;
							if (!(minEllipse[i].size.empty())) {
								if (minEllipse[i].size.area() > max_area) {

									max_area = minEllipse[i].size.area();
									best_contours[0] = contours[i];
									bestEllipse = minEllipse[i];
								}
							}
						}

						/// Draw contours + rotated rects + ellipses

						for (int i = 0; i < contours.size(); i++)
						{
							Scalar color = Scalar(224, 224, 224);
							// contour
							//drawContours(drawing, contours, i, color, 1, 8, vector<Vec4i>(), 0, Point());
							// ellipse
							ellipse(drawing, bestEllipse, color, 2, 8);
							// rotated rectangle
							//Point2f rect_points[4]; minRect[i].points(rect_points);
							//for (int j = 0; j < 4; j++)
								//line(drawing, rect_points[j], rect_points[(j + 1) % 4], color, 1, 8);
						}


						if (best_contours.size() > 0) {
							Scalar color = Scalar(224, 224, 224);

							//drawContours(drawing, best_contours, 0, color, 2, 8, hierarchy, 0, Point());

						}



						//imshow("contours", drawing);
						Mat gray_draw;
						cv::cvtColor(drawing, gray_draw, cv::COLOR_BGR2GRAY);

						Mat contrast;
						bitwise_or(gray, gray_draw, contrast);
						//imshow("contr", contrast);

						// Contours in ellipse
						vector<vector<Point> > contours3;
						Mat canny_output2;
						Mat contrast2;
						threshold(contrast, contrast2, 200, 255, cv::THRESH_BINARY);

						/// Detect edges using canny
						Canny(contrast2, canny_output2, 20, 20 * 2, 7);

						/// Find contours
						findContours(canny_output2, contours3, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, Point(0, 0));

						/// Draw contours

						drawing2 = Mat::zeros(canny_output2.size(), CV_8UC3);




						vector<vector<Point> >hull(contours3.size());
						for (size_t i = 0; i < contours3.size(); i++)
						{
							convexHull(contours3[i], hull[i]);
						}
						cout << contours3.size();
						for (size_t i = 0; i < contours3.size(); i++)
						{
							Scalar color = Scalar(192 + 10 * i, 192, 192);

							//drawContours(drawing2, contours3, (int)i, color);
							//drawContours(drawing2, hull, (int)1, color2,2,8,hierarchy,0,Point());

						}

						//imshow("contr", drawing2);

						//cout << n;
						vector<Point> half_hull0;

						for (int i = (hull[0].size()) / 2; i < hull[0].size(); i++) {

							half_hull0.push_back(hull[0][i]);

						}
						for (int i = 0 / 2; i < hull[0].size() / 2.; i++) {

							half_hull0.push_back(hull[0][i]);

						}
						for (int i = 0; i < half_hull0.size(); i = i + 1) {

							circle(drawing2, half_hull0[i], 2, Scalar(0, 0, 255), FILLED, LINE_AA);

							

						}
						circle(drawing2, Point(0, 0), 2, Scalar(0, 0, 255), FILLED, LINE_AA);
						int n = half_hull0.size();

						Mat final_draw;
						//cv::cvtColor(contrast, final_draw, cv::COLOR_GRAY2RGB);

						for (int i = 15; i < 130; i++) {

							for (int j = 41; j < 214; j++) {
								if (isInside(half_hull0, n - 1, Point(j, i)))//? cout << "Yes \n" : cout << "No \n";
								{
									contrast.at<uchar>(i, j) = 224;
								}

							}
						}
						//isInside(half_hull0, n-1, Point(0, 0)) ? cout << "Yes \n" : cout << "No \n";
						//imshow("contr", contrast);

						InImage_Ultra = contrast;



						//----------------------------------------------------------------------------------------------------//
						//---------------------------------------End of Image Segmentation-----------------------------------//
						//--------------------------------------------------------------------------------------------------//

					}

					//----- Put tumor in a larger image ------//

					Mat black = Mat::zeros(734, 1176, CV_8UC3);
					Mat another;
					Mat grayan;
					//cvtColor(InImage_Ultra, grayan, cv::COLOR_BGR2GRAY);
					grayan = InImage_Ultra;
					cvtColor(grayan, another, cv::COLOR_GRAY2RGB);
					Rect Rec3(1160, 285, 230, 170);

					Mat final_im = another;


					final_im.copyTo(black(cv::Rect(1160 - 614, 285 - 72, 230, 170)));

					imshow("Only Ultra", black);

					imwrite(name, black);
				}
				/// Show in a window
				//namedWindow("Contours", CV_WINDOW_AUTOSIZE);
				
				


				std::string tring = std::string(name);
				const char* fileName = tring.c_str();
				igsioStatus stat = igsioVideoFrame::ReadImageFromFile(video, fileName);

				trackedframe.SetImageData(video);


				///Matrices
				vtkMatrix4x4 *matrix4_Probe_Ref = vtkMatrix4x4::New();
				
				/*Markers[0].Tvec.ptr<float>(0)[0];
				Markers[0].Tvec.ptr<float>(0)[1];
				Markers[0].Tvec.ptr<float>(0)[2];

				Markers[0].Rvec.ptr<float>(0)[0];
				Markers[0].Rvec.ptr<float>(0)[1];
				Markers[0].Rvec.ptr<float>(0)[2];*/

				//matrix4_Probe_Ref->SetElement(0, 0, *R.ptr<double>(0, 0));
				
				for (int k = 0; k < 3; k++) {
					for (int j = 0; j < 3; j++) {
						matrix4_Probe_Ref->SetElement(k, j, *R.ptr<double>(k, j));
					}
				}
				matrix4_Probe_Ref->SetElement(0, 3, 1000.*tvecs[0][0]);
				matrix4_Probe_Ref->SetElement(1, 3, 1000.*tvecs[0][1]);
				matrix4_Probe_Ref->SetElement(2, 3, 1000.*tvecs[0][2]);
				matrix4_Probe_Ref->SetElement(3, 3, 1.);

				matrix4_Probe_Ref->SetElement(3, 0, 0.f);
				matrix4_Probe_Ref->SetElement(3, 1, 0.f);
				matrix4_Probe_Ref->SetElement(3, 2, 0.f);


				//1st column
				/*matrix4_Probe_Ref->SetElement(0, 0, Markers[0].Rvec.ptr<float>(0)[0]);
				matrix4_Probe_Ref->SetElement(1, 0, 0);
				matrix4_Probe_Ref->SetElement(2, 0, 0);
				matrix4_Probe_Ref->SetElement(3, 0, 0);

				//2nd column
				matrix4_Probe_Ref->SetElement(0, 1, 0);
				matrix4_Probe_Ref->SetElement(1, 1, Markers[0].Rvec.ptr<float>(0)[1]);
				matrix4_Probe_Ref->SetElement(2, 1, 0);
				matrix4_Probe_Ref->SetElement(3, 1, 0);

				//3rd column
				matrix4_Probe_Ref->SetElement(0, 2, 0);
				matrix4_Probe_Ref->SetElement(1, 2, 0);
				matrix4_Probe_Ref->SetElement(2, 2, Markers[0].Rvec.ptr<float>(0)[2]);
				matrix4_Probe_Ref->SetElement(3, 2, 0);

				//4th column
				matrix4_Probe_Ref->SetElement(0, 3, Markers[0].Tvec.ptr<float>(0)[0]);
				matrix4_Probe_Ref->SetElement(1, 3, Markers[0].Tvec.ptr<float>(0)[1]);
				matrix4_Probe_Ref->SetElement(2, 3, Markers[0].Tvec.ptr<float>(0)[2]);
				matrix4_Probe_Ref->SetElement(3, 3, 1);
				*/



				const igsioTransformName  nameProbRef = igsioTransformName(std::string("Probe"), std::string("Reference"));
				trackedframe.SetFrameTransform(nameProbRef, matrix4_Probe_Ref);
				trackedframe.SetTimestamp(timestamp);
				trackedframe.SetFrameTransformStatus(nameProbRef, ToolStatus::TOOL_OK);

				trackedFrameList->AddTrackedFrame(&trackedframe);


				timestamp = timestamp + 0.08;
				frames_counter++;

				detection_happened = false;
			}
			//cerr << "\nNumber of Correct data are:" << trackedFrameList->Size()<<endl;

			/*vtkMatrix4x4 *matrix4_Probe_Tracker = vtkMatrix4x4::New();
			vtkMatrix4x4 *matrix4_Ref_Tracker = vtkMatrix4x4::New();

			int elem = 0;
			for (int k = 0; k < 4; k++) {
				for (int j = 0; j < 4; j++) {
					matrix4_Probe_Tracker->SetElement(k, j, aray[i][elem]);
					elem++;
				}
			}
			elem = 0;

			for (int k = 0; k < 4; k++) {
				for (int j = 0; j < 4; j++) {
					matrix4_Ref_Tracker->SetElement(k, j, aray[i + 1][elem]);
					elem++;
				}
			}*/

			/*const igsioTransformName  nameProbTracker = igsioTransformName(std::string("Probe"), std::string("Tracker"));
			trackedframe.SetFrameTransform(nameProbTracker, matrix4_Probe_Tracker);
			const igsioTransformName  nameRefTracker = igsioTransformName(std::string("Reference"), std::string("Tracker"));
			trackedframe.SetFrameTransform(nameRefTracker, matrix4_Ref_Tracker);
			trackedframe.SetTimestamp(timestamp);
			trackedframe.SetFrameTransformStatus(nameProbTracker, ToolStatus::TOOL_OK);
			trackedframe.SetFrameTransformStatus(nameRefTracker, ToolStatus::TOOL_OK);
			//const igsioTransformName  nameRefTracker = igsioTransformName(std::string("Reference"));

			//trackedframe.SetFrameTransformStatus()
			trackedFrameList->AddTrackedFrame(&trackedframe);


			timestamp = timestamp + 0.08;*/

			/////////// End of frames ////////////////////////////////


		} while (frames_counter<used_frames && key != 27 && vreader.grab() && vreader_Ultra.grab());  // wait for esc to be pressed

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




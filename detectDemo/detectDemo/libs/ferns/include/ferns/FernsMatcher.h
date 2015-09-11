/*
* Filename: FernsMatcher.h
* Author:   Xiehong
* Date:     2013.2.7
*/

#include <stdio.h>
#include <stdlib.h>

#include "FernsDllExport.h"
#include "affine_transformation_range.h"


struct _IplImage;
typedef _IplImage IplImage;

class keypoint;
class affine_image_generator06;
class pyr_yape06;
class fine_gaussian_pyramid;
class fern_based_point_classifier;



// Parameters of ferns
struct FernsParams 
{
	affine_transformation_range range;
	int maximum_number_of_points_on_model;
	int number_of_generated_views;
	double minimum_number_of_views_rate;
	int patch_size;
	int yape_radius;
	int number_of_octaves;
	int number_of_ferns;
	int number_of_tests_per_fern;
	int number_of_samples_for_refinement;
	int number_of_samples_for_test;
	int roi_u[4];	// ul, ur, br, bl
	int roi_v[4];

	// Default data
	FernsParams()
	{
		maximum_number_of_points_on_model = 400;
		number_of_generated_views = 5000;
		minimum_number_of_views_rate = 0.0;
		patch_size = 32;
		yape_radius = 7;
		number_of_octaves = 4;
		number_of_ferns = 30;
		number_of_tests_per_fern = 12;
		number_of_samples_for_refinement = 10000;
		number_of_samples_for_test = 200;
		for (int i=0; i<4; i++)
		{
			roi_u[i] = -1;
			roi_v[i] = -1;
		}
	}
};



// Class of ferns training and matching
class FERNS_CLASS FernsMatcher
{
public:
	FernsMatcher();
	virtual ~FernsMatcher();

public:
	struct Pointf
	{
		float u;
		float v;

		Pointf(float _u = 0, float _v = 0) : u(_u), v(_v) {}
	};

public:
	// Train ferns based classifier
	bool Train(IplImage *pGrayModel, const FernsParams &params, const char *modelPointsFilename = NULL);

	// Load trained classifier
	bool Load(const char *filename);

	// Save trained classifier
	bool Save(const char *filename);

	// Set the maximum number of points we want to detect
	void SetMaxNumberOfPointsDetect(int nMax);

	// Detect features
	int Detect(const IplImage *pGrayImage);

	// Match features
	int Match(float minClassScore = 0);

	// Get matched points
	void GetMatches(Pointf *p1, Pointf *p2, float *scores = NULL);

	void GetModelROI(int u[4], int v[4]);

	int GetModelPointsCount() { return number_of_model_points; }
	int GetDetectPointsCount() { return number_of_detected_points; }
	int GetMatchesCount() { return number_of_matches; }

protected:
	void GetBoundRect(int u[4], int v[4], int &x_min, int &y_min, int &x_max, int &y_max);
	void SaveImageOfModelPoints(const char *filename, int patch_size);

	void DetectMostStableModelPoints(int maximum_number_of_points_on_model, int yape_radius, int number_of_octaves,
									int number_of_generated_images, double minimum_number_of_views_rate);

protected:
	fern_based_point_classifier *classifier;
	affine_image_generator06 *image_generator;
	pyr_yape06 *point_detector;
	fine_gaussian_pyramid *pyramid;

	IplImage *model_image;

	keypoint *model_points;
	keypoint *detected_points;
	int number_of_model_points;
	int number_of_detected_points;
	int maximum_number_of_points_to_detect;

	int u_corner[4], v_corner[4];
	int number_of_matches;

	int patch_size;
	int yape_radius;
	int number_of_octaves;

	float mean_recognition_rate;
};
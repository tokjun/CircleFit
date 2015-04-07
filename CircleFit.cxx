/*=========================================================================

  Program:   Fit 2D plane to a set of 3D points
  Language:  C++
  Author:    Junichi Tokuda, Ph.D. (Brigham and Women's Hospital)

  Copyright (c) Brigham and Women's Hospital. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/


#include <complex>
#include "itkListSample.h"
#include "itkCovarianceSampleFilter.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkAffineTransform.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

typedef std::vector<double> CoordType;
typedef std::vector<CoordType> CoordSetType;

typedef itk::Vector< double, 3 > VectorType;
typedef itk::Point< double, 3> PointType;
typedef itk::Matrix< double, 3, 3 > MatrixType;
typedef itk::FixedArray< double, 3 > ArrayType;
typedef itk::AffineTransform< double, 3 > TransformType;
typedef itk::Statistics::ListSample< VectorType > PointListType;
typedef PointListType::Iterator PointListIteratorType;
typedef itk::Statistics::CovarianceSampleFilter< PointListType > CovarianceAlgorithmType;
typedef itk::SymmetricEigenAnalysis< CovarianceAlgorithmType::MatrixType, ArrayType, MatrixType > SymmetricEigenAnalysisType;

int CalcIntersectionOfPerpendicularBisectors2D(VectorType& p1, VectorType& p2, VectorType& p3, VectorType& intersec, double radius);
double CalcAverageMinDistanceOfRotatedPoints(PointListType::Pointer rotatingPoints, PointListType::Pointer fixedPoints,
					     VectorType& principalVector, int angle,
					     MatrixType& rotationMatrix);
  
// Load points from the csv file. If successful, return > 0
int LoadPoints(const char* filename, PointListType* points)
{
  std::string line;
  std::ifstream configfile (filename);
  if (!configfile.is_open())
    {
    std::cerr << "Unable to open file" << std::endl; 
    return 0;
    }

  points->Clear();

  VectorType point;

  std::cout << "Marker Configuration File: " << filename << std::endl;
  while ( configfile.good() )
    {
    int pointId = 0;
    while(std::getline(configfile,line))
      {
      std::vector<double> values;
      std::stringstream  lineStream(line);
      std::string        cell;
      values.clear();
      while(std::getline(lineStream,cell,','))
        {
        values.push_back(::atof(cell.c_str()));
        }
      if (values.size() == 3)
        {
        point[0] = values[0];
        point[1] = values[1];
        point[2] = values[2];
        points->PushBack(point);
        std::cout << "Fiducial #" << pointId << ": "
                  << "Point=("  << point[0] << ", " << point[1] << ", " << point[2] << "); "
                  << std::endl;
        pointId ++;
        } 
      else
        {
        if (values.size() != 0)
          {
          std::cout << "Invalid format!" << std::endl;
          return 0;
          }
        }
      }
    }
  std::cout << std::endl;
  configfile.close();
  
  return 1;
}


int main( int argc, char * argv [] )
{

  if( argc < 2 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " srcPoints dstPoints [radius]" << std::endl;
    return EXIT_FAILURE;
    }
  const unsigned int Dimension = 2;

  double radius = -1.0;
  if (argc > 3)
    {
    radius = atof(argv[3]);
    }

  //----------------------------------------
  // Load the points

  PointListType::Pointer srcPoints, dstPoints;
  srcPoints = PointListType::New();
  dstPoints = PointListType::New();

  LoadPoints(argv[1], srcPoints);
  LoadPoints(argv[2], dstPoints);

  //----------------------------------------
  // Perform PCA

  CovarianceAlgorithmType::Pointer covarianceAlgorithm = 
    CovarianceAlgorithmType::New();
  
  covarianceAlgorithm->SetInput( dstPoints );
  covarianceAlgorithm->Update();
    
  // Perform Symmetric Eigen Analysis
  SymmetricEigenAnalysisType analysis ( 3 );
  ArrayType eigenValues;
  MatrixType eigenMatrix;
  analysis.SetOrderEigenMagnitudes( true );
  analysis.ComputeEigenValuesAndVectors( covarianceAlgorithm->GetCovarianceMatrix(),
                                         eigenValues, eigenMatrix );    

  //----------------------------------------
  // Extract the normal vector 
  // The first eigen vector (with the minimal eigenvalue) is
  // the normal vector of the fitted plane.

  VectorType nx =  eigenMatrix[2];
  VectorType ny =  eigenMatrix[1];
  VectorType nz =  eigenMatrix[0];
  VectorType principalVector = nz;

  //----------------------------------------
  // Calculate the average position of the all points.
  // This is used as the origin of the new coordinate system on the
  // fitted plane.

  CovarianceAlgorithmType::MeasurementVectorType meanPoint;
  meanPoint = covarianceAlgorithm->GetMean();
  
  //----------------------------------------
  // Project all the points to the fitted plane.

  PointListType::Pointer projectedPoints = PointListType::New();
  for (PointListIteratorType iter = dstPoints->Begin(); iter != dstPoints->End(); ++iter)
    {
    VectorType p1;
    VectorType p2;
    p1 = iter.GetMeasurementVector() - meanPoint;
    p2[0] = p1*nx;
    p2[1] = p1*ny;
    p2[2] = 0.0;
    projectedPoints->PushBack(p2);
    }

  //----------------------------------------
  // Pick up every combination of three points from the list and calculate
  // the intersection of the perpendicular bisectors of the two chords connecting
  // the three points.

  VectorType meanIntersect;
  meanIntersect[0] = 0.0;
  meanIntersect[1] = 0.0;
  meanIntersect[2] = 0.0;

  int nPoints = 0;
  int nPointsUsed = 0;

  for (PointListIteratorType iter1 = projectedPoints->Begin(); iter1 != projectedPoints->End(); ++iter1)
    {
    PointListIteratorType iter2 = iter1;
    for (++iter2; iter2 != projectedPoints->End(); ++iter2)
      {
      PointListIteratorType iter3 = iter2;
      for (++iter3; iter3 != projectedPoints->End(); ++iter3)
        {
        VectorType p1 = iter1.GetMeasurementVector();
        VectorType p2 = iter2.GetMeasurementVector();
        VectorType p3 = iter3.GetMeasurementVector();
        VectorType c;

        nPoints ++;
        if (CalcIntersectionOfPerpendicularBisectors2D(p1, p2, p3, c, radius) > 0)
          {
          meanIntersect = meanIntersect + c;
          nPointsUsed ++;
          }
        }
      }
    }
  meanIntersect = meanIntersect / nPointsUsed;

  //----------------------------------------
  // Transform the center point to the original coordinate system

  VectorType center = meanIntersect[0] * nx + meanIntersect[1] * ny + meanIntersect[2] * nz + meanPoint;

  std::cout << "Number of estimated center points: " << nPoints << std::endl;
  std::cout << "Number of estimated center points used: " << nPointsUsed << std::endl;
  std::cout << "Center = " << center << std::endl;  
  std::cout << std::endl;

  //----------------------------------------
  // Calculate matrix from original coordinate system to plane coordinate system

  MatrixType originalToPlaneMatrix;
  for (int i = 0; i < 3; ++i)
    {
    originalToPlaneMatrix[i][0] = nx[i];
    originalToPlaneMatrix[i][1] = ny[i];
    originalToPlaneMatrix[i][2] = nz[i];
    }
    
  //----------------------------------------
  // Rotate points from original position to in-plane position
  
  PointListType::Pointer inPlanePoints = PointListType::New();
  for (PointListIteratorType iter = srcPoints->Begin(); iter != srcPoints->End(); ++iter)
    {
    VectorType pp = originalToPlaneMatrix*iter.GetMeasurementVector() + center;
    inPlanePoints->PushBack(pp);
    }

  //----------------------------------------
  // Rotate point around principal vector and compute average minimum distance
  
  double globalAverageMinimumDistance = -1.0;
  MatrixType tmpRotationMatrix;
  MatrixType globalRotationMatrix;
  double globalRotationAngle = 0.0;
  bool circleFlipped = false;

  // TODO: Optimize. Use a step of 3 degrees (for example) and find minimum average distance.
  //       When found, use a step of 1 degree between the 3 degrees previously found.
  //       Step of 1 = 360 iterations
  //       Optimized (step of 3) = (360 / 3) + 3 = 123 iterations
  //       Optimized (step of 10) = (360 / 10) + 10 = 46 iterations
  //       If step is too big though, average minimum can be missed, and registration fail
  for (int angle = 0; angle < 360; ++angle)
    {
    double averageMinDistance = CalcAverageMinDistanceOfRotatedPoints(inPlanePoints, dstPoints,
								      principalVector, angle, 
								      tmpRotationMatrix);

    if (globalAverageMinimumDistance < 0 || averageMinDistance < globalAverageMinimumDistance)
      {
      globalAverageMinimumDistance = averageMinDistance;
      globalRotationMatrix = tmpRotationMatrix;
      globalRotationAngle = angle;
      }
    }

  //----------------------------------------
  // Rotate circle around one of the other axis, nx or ny, and recalculate average minimum distance
  // to also find the global minimum (including symmetry)

  TransformType::Pointer flippingTransform = TransformType::New();
  flippingTransform->Rotate3D(nx, M_PI);

  PointListType::Pointer flippedInPlanePoints = PointListType::New();
  for (PointListIteratorType iter = inPlanePoints->Begin(); iter != inPlanePoints->End(); ++iter)
    {
    flippedInPlanePoints->PushBack(flippingTransform->GetMatrix()*iter.GetMeasurementVector());
    }

  for (int angle = 0; angle < 360; ++angle)
    {
    double averageMinDistance = CalcAverageMinDistanceOfRotatedPoints(flippedInPlanePoints, dstPoints,
								      principalVector, angle,
								      tmpRotationMatrix);

    if (globalAverageMinimumDistance < 0 || averageMinDistance < globalAverageMinimumDistance)
      {
      globalAverageMinimumDistance = averageMinDistance;
      globalRotationMatrix = tmpRotationMatrix;
      globalRotationAngle = angle;
      circleFlipped = true;
      }
    }

  std::cout << "Best Fitting angle: " << globalRotationAngle << " (" << (circleFlipped ? "Flipped" : "Not flipped") << ")" << std::endl
	    << "Average Minimum Distance: " << globalAverageMinimumDistance << std::endl;

  MatrixType registrationMatrix = circleFlipped ? 
    globalRotationMatrix*flippingTransform->GetMatrix()*originalToPlaneMatrix :
    globalRotationMatrix*originalToPlaneMatrix;
  std::cout << "Registration Matrix: " << std::endl << registrationMatrix << std::endl;

  //----------------------------------------
  // Output registered points for checking
  

  std::ofstream outputFile("registeredPoints.csv");
  for (PointListIteratorType iter = srcPoints->Begin(); iter != srcPoints->End(); ++iter)
    {
    VectorType tmpVector = center + registrationMatrix*iter.GetMeasurementVector();
    outputFile << tmpVector[0] << "," << tmpVector[1] << "," << tmpVector[2] << std::endl;
    }
  outputFile.close();

  return EXIT_SUCCESS;
}

//--------------------------------------------------------------------------------
// Return 0, if the data set may not give a good estimate.
int CalcIntersectionOfPerpendicularBisectors2D(VectorType& p1, VectorType& p2, VectorType& p3, VectorType& intersect, double radius)
{

  const double posErr = 5.0;

  // Compute the bisecting points between p1 and p2 (m1), and p2 and p3 (m2)
  VectorType m1 = (p1+p2)/2.0;
  VectorType m2 = (p2+p3)/2.0;

  // Compute the normal vectors along the perpendicular bisectors
  VectorType v12 = (p2-p1);
  VectorType v23 = (p3-p2);

  if (v12.GetNorm() < posErr || v23.GetNorm() < posErr)
    {
    return 0;
    }

  v12.Normalize();
  v23.Normalize();

  VectorType n1;
  n1[0] = -v12[1];
  n1[1] = v12[0];
  n1[2] = v12[2];

  VectorType n2;
  n2[0] = -v23[1];
  n2[1] = v23[0];
  n2[2] = v23[2];
  
  // Compute the projection of m2 onto the perpendicular bisector of p1p2
  VectorType h = m1 + ((m2-m1)*n1)*n1;

  // The intersecting point of the two perpendicular bisectors (= estimated
  // center of fitted circle) is 'c' can be written as:
  //
  //    <c> = <m2> + a * <n2>
  //
  // where 'a' is a scalar value. Projection of 'c' on the m2h is 'h'
  //
  //    a * <n2> * (<h> - <m2>)/(|<h> - <m2>|) = |<h> - <m2>|
  //    a = |<h> - <m2>|^2 / {<n2> * (<h> - <m2>)}
  //

  VectorType m2h = h-m2;
  VectorType::RealValueType a = m2h.GetSquaredNorm() / (n2 * m2h);
  
  intersect = m2 + a * n2;

  // Validate by radius, if radius is greater than 0
  if (radius > 0)
    {
    VectorType d1 = p1-intersect;
    VectorType d2 = p2-intersect;
    VectorType d3 = p3-intersect;
    
    if (fabs(d1.GetNorm()-radius) > posErr ||
        fabs(d2.GetNorm()-radius) > posErr ||
        fabs(d3.GetNorm()-radius) > posErr)
      {
      return 0;
      }
    }

  return 1;
}

//--------------------------------------------------------------------------------
// Return average minimum distance between rotated and fixed pointsets,
// as well as the rotation matrix used
double CalcAverageMinDistanceOfRotatedPoints(PointListType::Pointer rotatingPoints, PointListType::Pointer fixedPoints,
					   VectorType& principalVector, int angle, 
					   MatrixType& rotationMatrix)
{
  // Rotate points
  TransformType::Pointer rotationTransform = TransformType::New();
  rotationTransform->Rotate3D(principalVector, angle * M_PI / 180);
  rotationMatrix = rotationTransform->GetMatrix();
  
  double averageMinDistance = 0.0;
  int numberOfPoints = 0;
  for (PointListIteratorType iter1 = rotatingPoints->Begin(); iter1 != rotatingPoints->End(); ++iter1)
    {
    // Rotate point
    VectorType rp = rotationMatrix*iter1.GetMeasurementVector();

    // Compute minDistance
    double minDistance = -1.0;
    for (PointListIteratorType iter2 = fixedPoints->Begin(); iter2 != fixedPoints->End(); ++iter2)
      {
      VectorType fp = iter2.GetMeasurementVector();
      
      double distance = std::sqrt(std::pow(rp[0]-fp[0],2) + 
				  std::pow(rp[1]-fp[1],2) + 
				  std::pow(rp[2]-fp[2],2));
      if (minDistance < 0 || distance < minDistance)
	{
	minDistance = distance;
	}
      }
    averageMinDistance += minDistance;
    numberOfPoints++;
    }
  averageMinDistance /= numberOfPoints;
  return averageMinDistance;
}

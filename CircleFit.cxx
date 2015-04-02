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
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

typedef std::vector<double> CoordType;
typedef std::vector<CoordType> CoordSetType;

typedef itk::Vector< double, 3 > VectorType;
typedef itk::Statistics::ListSample< VectorType > PointListType;


int CalcIntersectionOfPerpendicularBisectors2D(VectorType& p1, VectorType& p2, VectorType& p3, VectorType& intersec, double radius);
  
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
  configfile.close();
  
  return 1;
}


int main( int argc, char * argv [] )
{

  if( argc < 2 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " pointData [radius]" << std::endl;
    return EXIT_FAILURE;
    }
  const unsigned int Dimension = 2;

  typedef std::complex< float >              PixelType;

  PointListType::Pointer points;
  points = PointListType::New();

  double radius = -1.0;
  if (argc > 2)
    {
    radius = atof(argv[2]);
    }

  //----------------------------------------
  // Load the points
  LoadPoints(argv[1], points);

  //----------------------------------------
  // Perform PCA
  typedef itk::Statistics::CovarianceSampleFilter< PointListType > 
    CovarianceAlgorithmType;
  CovarianceAlgorithmType::Pointer covarianceAlgorithm = 
    CovarianceAlgorithmType::New();
  
  covarianceAlgorithm->SetInput( points );
  covarianceAlgorithm->Update();
  
  //std::cout << "Sample covariance = " << std::endl ; 
  //std::cout << covarianceAlgorithm->GetCovarianceMatrix() << std::endl;
  
  // Perform Symmetric Eigen Analysis
  typedef itk::FixedArray< double, 3 > EigenValuesArrayType;
  typedef itk::Matrix< double, 3, 3 > EigenVectorMatrixType;
  typedef itk::SymmetricEigenAnalysis< CovarianceAlgorithmType::MatrixType,
    EigenValuesArrayType, EigenVectorMatrixType > SymmetricEigenAnalysisType;
  SymmetricEigenAnalysisType analysis ( 3 );
  EigenValuesArrayType eigenValues;
  EigenVectorMatrixType eigenMatrix;
  analysis.SetOrderEigenMagnitudes( true );
  analysis.ComputeEigenValuesAndVectors( covarianceAlgorithm->GetCovarianceMatrix(),
                                         eigenValues, eigenMatrix );    

  // Print out the result of PCA
  std::cout << eigenValues << std::endl;
  std::cout << eigenMatrix << std::endl;


  //----------------------------------------
  // Extract the normal vector 
  // The first eigen vector (with the minimal eigenvalue) is
  // the normal vector of the fitted plane.

  VectorType principalVector = eigenMatrix[0];
  std::cout << principalVector << std::endl;

  //----------------------------------------
  // Calculate the average position of the all points.
  // This is used as the origin of the new coordinate system on the
  // fitted plane.

  CovarianceAlgorithmType::MeasurementVectorType meanPoint;
  meanPoint = covarianceAlgorithm->GetMean();
  std::cout << "Sample mean = " << meanPoint << std::endl ; 
  
  //----------------------------------------
  // Project all the points to the fitted plane.

  PointListType::Pointer transPoints = PointListType::New();
  typedef PointListType::Iterator IteratorType;
  IteratorType iter = points->Begin();

  VectorType nx =  eigenMatrix[2];
  VectorType ny =  eigenMatrix[1];
  VectorType nz =  eigenMatrix[0];
  
  while (iter != points->End())
    {
    VectorType p1;
    VectorType p2;
    p1 = iter.GetMeasurementVector() - meanPoint;
    p2[0] = p1*nx;
    p2[1] = p1*ny;
    //p2[2] = p1*nz;
    p2[2] = 0.0;
    transPoints->PushBack(p2);
    ++ iter;
    }

  // Pick up every combination of three points from the list and calculate
  // the intersection of the perpendicular bisectors of the two chords connecting
  // the three points.

  VectorType meanIntersect;
  meanIntersect[0] = 0.0;
  meanIntersect[1] = 0.0;
  meanIntersect[2] = 0.0;

  int nPoints = 0;
  int nPointsUsed = 0;

  for (IteratorType iter1 = transPoints->Begin(); iter1 != transPoints->End(); ++ iter1)
    {
    IteratorType iter2 = iter1;
    for (++ iter2; iter2 != transPoints->End(); ++ iter2)
      {
      IteratorType iter3 = iter2;
      for (++ iter3; iter3 != transPoints->End(); ++ iter3)
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

  // Transform the center point to the original coordinate system
  VectorType center = meanIntersect[0] * nx + meanIntersect[1] * ny + meanIntersect[2] * nz + meanPoint;

  std::cout << "Number of points: " << nPoints << std::endl;
  std::cout << "Number of points used: " << nPointsUsed << std::endl;
  std::cout << "Center = " << center << std::endl;

  return EXIT_SUCCESS;
}


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

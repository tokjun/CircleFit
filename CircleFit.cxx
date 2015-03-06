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
    std::cerr << argv[0] << " pointData " << std::endl;
    return EXIT_FAILURE;
    }
  const unsigned int Dimension = 2;

  typedef std::complex< float >              PixelType;

  PointListType::Pointer points;
  points = PointListType::New();

  LoadPoints(argv[1], points);

  typedef itk::Statistics::CovarianceSampleFilter< PointListType > 
    CovarianceAlgorithmType;
  CovarianceAlgorithmType::Pointer covarianceAlgorithm = 
    CovarianceAlgorithmType::New();
  
  covarianceAlgorithm->SetInput( points );
  covarianceAlgorithm->Update();
  
  //std::cout << "Sample covariance = " << std::endl ; 
  //std::cout << covarianceAlgorithm->GetCovarianceMatrix() << std::endl;
  
  CovarianceAlgorithmType::MeasurementVectorType meanVector;
  meanVector = covarianceAlgorithm->GetMean();
  std::cout << "Sample mean = " << meanVector << std::endl ; 
  
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

  std::cerr << eigenValues << std::endl;
  std::cerr << eigenMatrix << std::endl;

  // The third eigen vector is the normal vector of the fitted plane
  VectorType principalVector = eigenMatrix[0];
  std::cerr << principalVector << std::endl;

  // Transform each point

  PointListType::Pointer transPoints = PointListType::New();

  typedef PointListType::Iterator IteratorType;
  IteratorType iter = points->Begin();

  while (iter != points->End())
    {
    VectorType p1;
    VectorType p2;

    p1 = iter.GetMeasurementVector() - meanVector;

    p2[0] = p1*eigenMatrix[2];
    p2[1] = p1*eigenMatrix[1];
    p2[2] = p1*eigenMatrix[0];
    transPoints->PushBack(p2);

    std::cout << "P2=("  << p2[0] << ", " << p2[1] << ", " << p2[2] << "); " << std::endl;

    ++ iter;
    }

  // Least square fit
  


  return EXIT_SUCCESS;
}

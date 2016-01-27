/*=========================================================================

  Program:   3D Slicer Marker Registration CLI
  Module:    LineMarkerRegistration
  Language:  C++
  Author:    Junichi Tokuda, Ph.D. (Brigham and Women's Hospital)

  Copyright (c) Brigham and Women's Hospital. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <algorithm>

#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkLabelToLineImageFilter.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkTransformFileWriter.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkEuclideanDistanceLineMetric.h"
#include "itkChangeLabelImageFilter.h"

#include "itkAffineTransform.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

#include "itkEuler3DTransform.h"
#include "itkRigid3DTransform.h"
#include "itkTranslationTransform.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "itkPointSetToPointSetRegistrationMethod.h"

#include "itkPluginUtilities.h"
#include "LineMarkerRegistrationCLP.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace {

class CommandIterationUpdate : public itk::Command 
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {};

public:

  typedef itk::LevenbergMarquardtOptimizer     OptimizerType;
  typedef const OptimizerType *                OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    OptimizerPointer optimizer = 
                         dynamic_cast< OptimizerPointer >( object );

    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }

    }
   
};


typedef std::vector<double> CoordType;
typedef std::vector<CoordType> CoordSetType;


// Load the configuration file. If successful, return > 0
// TODO: This is a tentative solution. The fucntion has to be replaced by XML parser.
int LoadMarkerConfiguration(const char* filename, CoordSetType& points)
{
  std::string line;
  std::ifstream configfile (filename);
  if (!configfile.is_open())
    {
    std::cerr << "Unable to open file" << std::endl; 
    return 0;
    }

  points.clear();

  CoordType point(3);
  CoordType norm(3);

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
      if (values.size() == 6)
        {
        point[0] = values[0];
        point[1] = values[1];
        point[2] = values[2];
        norm[0] = values[3];
        norm[1] = values[4];
        norm[2] = values[5];
        points.push_back(point);
        points.push_back(norm);
        std::cout << "Fiducial #" << pointId << ": "
                  << "Point=("  << point[0] << ", " << point[1] << ", " << point[2] << "); "
                  << "Normal=(" << norm[0] << ", " << norm[1] << ", " << norm[2] << ")"
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


// For sorting
typedef struct {
  int label;
  double meanVesselness;
} LabelVesselness;

bool SortByVesselness (const LabelVesselness& lhs, const LabelVesselness &rhs) { return lhs.meanVesselness > rhs.meanVesselness; }
  
  
template<class T> int DoIt( int argc, char * argv[], T )
{
  PARSE_ARGS;

  const     unsigned int        Dimension       = 3;

  typedef   T                   FileInputPixelType;
  typedef   float               InternalPixelType;
  typedef   int                 OutputPixelType;


  typedef   itk::Image< FileInputPixelType, Dimension > FileInputImageType;
  typedef   itk::Image< InternalPixelType, Dimension >  InternalImageType;
  typedef   itk::Image< OutputPixelType, Dimension >    OutputImageType;

  typedef   itk::ImageFileReader< InternalImageType >  ReaderType;
  typedef   itk::ImageFileWriter< OutputImageType > WriterType;
  //typedef   itk::ImageFileWriter< InternalImageType > WriterType;

  // Smoothing filter
  typedef   itk::SmoothingRecursiveGaussianImageFilter<
  InternalImageType, InternalImageType > SmoothingFilterType;
  
  // Line enhancement filter
  typedef   itk::Hessian3DToVesselnessMeasureImageFilter<
    InternalPixelType > LineFilterType;
  //typedef   itk::Hessian3DToNeedleImageFilter<
  //  InternalPixelType > LineFilterType;

  typedef   itk::ConnectedComponentImageFilter<
    OutputImageType, OutputImageType >  CCFilterType;
  typedef   itk::RelabelComponentImageFilter<
    OutputImageType, OutputImageType > RelabelType;

  // Line detection filter
  typedef   itk::LabelToLineImageFilter<
  OutputImageType, OutputImageType > LabelLineFilterType;

  // Declare the type of enhancement filter - use ITK's 3D vesselness (Sato)
  //typedef itk::Hessian3DToVesselnessMeasureImageFilter<double> VesselnessFilterType;
  // Declare the type of multiscale enhancement filter
  typedef itk::RescaleIntensityImageFilter<InternalImageType> RescaleFilterType;
  typedef itk::Hessian3DToVesselnessMeasureImageFilter<double> VesselnessFilterType;

  typedef itk::NumericTraits< InternalPixelType >::RealType RealPixelType;
  typedef itk::SymmetricSecondRankTensor< RealPixelType, Dimension > HessianPixelType;
  typedef itk::Image< HessianPixelType, Dimension >                  HessianImageType;
  typedef itk::MultiScaleHessianBasedMeasureImageFilter< InternalImageType, HessianImageType, InternalImageType >
    MultiScaleEnhancementFilterType;

  typedef itk::BinaryThresholdImageFilter <InternalImageType, OutputImageType>
    BinaryThresholdImageFilterType;
 
  typename ReaderType::Pointer reader = ReaderType::New();  
  typename WriterType::Pointer writer = WriterType::New();
  typename SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();

  typename LineFilterType::Pointer lineFilter = LineFilterType::New();
  typename CCFilterType::Pointer CCFilter = CCFilterType::New();
  typename RelabelType::Pointer RelabelFilter = RelabelType::New();
  typename LabelLineFilterType::Pointer labelLineFilter = LabelLineFilterType::New();

  // Detect lines from the label map
  typedef itk::PointSet< float, 3 >   PointSetType;
  typedef typename LabelLineFilterType::LineTransformType TransformType;
  typedef itk::EuclideanDistanceLineMetric<PointSetType, PointSetType> LineDistanceMetric;
  
  typedef typename LineDistanceMetric::MovingPointSetType MovingPointSetType;
  typedef typename MovingPointSetType::PointsContainer MovingPointSetContainer;
  typedef typename MovingPointSetType::PointType MovingPointType;
  
  typedef typename LineDistanceMetric::FixedPointSetType FixedPointSetType;
  typedef typename FixedPointSetType::PointsContainer FixedPointSetContainer;
  typedef typename FixedPointSetType::PointType FixedPointType;

  typedef itk::LabelStatisticsImageFilter< InternalImageType, OutputImageType > LabelStatisticsType;

  //
  // Create pointset of the Z-frame
  // TODO: this has to be configurable (e.g. XML)
  //
  MovingPointSetType::Pointer movingPointSet = MovingPointSetType::New();
  MovingPointSetContainer::Pointer movingPointSetContainer = MovingPointSetContainer::New();

  CoordSetType points;
  if (LoadMarkerConfiguration(markerConfigFile.c_str(), points) <= 0)
    {
    return EXIT_FAILURE ;
    }

  int id = 0;
  CoordSetType::iterator iter;
  for (iter = points.begin(); iter != points.end(); iter ++)
    {
    MovingPointType point;
    point[0] = (*iter)[0];
    point[1] = (*iter)[1];
    point[2] = (*iter)[2];
    movingPointSetContainer->InsertElement(id, point);
    id ++;
    }

  // Number of line markers in the frame. ('point' contains array of [point0, norm0, point1, norm1, ...]
  int nMarkers = points.size() / 2;
  
  reader->SetFileName( inputVolume.c_str() );
  writer->SetFileName( outputVolume.c_str() );

  smoothing->SetInput( reader->GetOutput() );
  smoothing->SetSigma( static_cast< double >(sigma1) );

  lineFilter->SetAlpha1( static_cast< double >(alpha1));
  lineFilter->SetAlpha2( static_cast< double >(alpha2));

  // We use only a signle-scale filter.
  MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter = MultiScaleEnhancementFilterType::New();
  multiScaleEnhancementFilter->SetInput(smoothing->GetOutput());
  multiScaleEnhancementFilter->SetSigmaMinimum(sigma2);
  multiScaleEnhancementFilter->SetSigmaMaximum(sigma2);
  multiScaleEnhancementFilter->SetNumberOfSigmaSteps(1);
  multiScaleEnhancementFilter->SetHessianToMeasureFilter (lineFilter);

  BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();
  thresholdFilter->SetInput(multiScaleEnhancementFilter->GetOutput());
  thresholdFilter->SetLowerThreshold(lowerThreshold);
  thresholdFilter->SetUpperThreshold(upperThreshold);
  thresholdFilter->SetInsideValue(255);
  thresholdFilter->SetOutsideValue(0);

  CCFilter->SetInput (thresholdFilter->GetOutput());
  CCFilter->FullyConnectedOff();

  try
    {
    CCFilter->Update();
    }
  catch (itk::ExceptionObject &err)
    {
    std::cerr << err << std::endl;
    return EXIT_FAILURE ;
    }

  RelabelFilter->SetInput ( CCFilter->GetOutput() );

  // Calculate minimum number of voxels in each object
  OutputImageType::SpacingType spacing = CCFilter->GetOutput()->GetSpacing();
  int min = (int) (minimumObjectSize/(spacing[0]*spacing[1]*spacing[2]));
  RelabelFilter->SetMinimumObjectSize( min );
  RelabelFilter->Update();

  MovingPointSetContainer::Pointer movingPointContainer = MovingPointSetContainer::New();

  labelLineFilter->SetInput( RelabelFilter->GetOutput() );
  //typedef typename RelabelType::ObjectSizeInPhysicalUnitsContainerType PhysicalSizeContainerType; // for ITK v4
  typedef std::vector<float> PhysicalSizeContainerType;
  PhysicalSizeContainerType objectSize = RelabelFilter->GetSizeOfObjectsInPhysicalUnits();
  int nObjects = RelabelFilter->GetNumberOfObjects();

  // Setup a map for change label filter
  typedef itk::ChangeLabelImageFilter<OutputImageType, OutputImageType> ChangeLabelFilter;
  typedef typename ChangeLabelFilter::ChangeMapType ChangeMapType;

  ChangeMapType changeMap;
  FixedPointType Centroid;

  Centroid[0] = 0.0;
  Centroid[1] = 0.0;
  Centroid[2] = 0.0;

  FixedPointSetType::Pointer fixedPointSet = FixedPointSetType::New();
  FixedPointSetContainer::Pointer fixedPointSetContainer = FixedPointSetContainer::New();
  unsigned int pointId = 0;

  LabelStatisticsType::Pointer labelStatistics = LabelStatisticsType::New();
  labelStatistics->SetLabelInput( RelabelFilter->GetOutput() );
  labelStatistics->SetInput( multiScaleEnhancementFilter->GetOutput() );
  labelStatistics->Update();
  
  std::cout << "Number of labels: " << labelStatistics->GetNumberOfLabels() << std::endl;
  std::cout << "Filter by Dimensions: " << FilterByDimensions << std::endl;
  std::cout << "Filter by Vesselness: " << FilterByVesselness << std::endl;

  // Sort the labels based on the mean 'vesselness' values
  std::vector< LabelVesselness > labels;
  for (int i = 0; i < nObjects; i ++)
    {
    LabelVesselness lv;
    lv.label = i+1;
    lv.meanVesselness = labelStatistics->GetMean(lv.label);
    labels.push_back(lv);
    }

  std::sort(labels.begin(), labels.end(), SortByVesselness);
  
  // Validate the detected objects in the order of mean vesselness value based on their diemensions.
  // if 'FilterByVesselness' is specified, the validation process will stop as soon as it finds N
  // valid objects, where N is the number of markers in the fiducial frame model.
  
  int nValidMarkers = 0;
  
  std::vector< LabelVesselness >::iterator lviter;
  for (lviter = labels.begin(); lviter != labels.end(); lviter ++)
    {
    bool fExclude = false;
      
    // According to ITK's manual:
    // "Once all the objects are relabeled, the application can query the number of objects and
    //  the size of each object. Object sizes are returned in a vector. The size of the background
    //  is not calculated. So the size of object #1 is GetSizeOfObjectsInPixels()[0], the size of object
    //  #2 is GetSizeOfObjectsInPixels()[1], etc."
    //float size = objectSize[i];
    //int label = i + 1;
    int label = lviter->label;
    float size = objectSize[label-1];

    // NOTE: size < minimumObjectSize might be redundant here, since it was already applied in the RelabelFilter.
    if (size < minimumObjectSize || size > maximumObjectSize)
      {
      // Out of size criteria
      fExclude = true;
      }
    else
      {
      // TODO: It is probably a good idea to also check the length and thickness of
      // the segmented area using LabelToLineImageFilter::GetAxisLength().
      FixedPointType point;
      FixedPointType norm;
      labelLineFilter->SetLabel( label );
      labelLineFilter->Update();
      typedef typename LabelLineFilterType::VectorType VectorType;
      VectorType axisLength;
      labelLineFilter->GetAxisLength(axisLength);
      
      // If the principal axis is less than the minimumPrincipalAxisLength or 
      // if any of the minor axes are longer than maximumMinorAxis
      if (FilterByDimensions &&
          (axisLength[0] < minimumPrincipalAxisLength ||
           axisLength[0] > maximumPrincipalAxisLength ||
           axisLength[1] > maximumMinorAxis ||
           axisLength[2] > maximumMinorAxis))
        {
        fExclude = true;
        }
      else
        {
        nValidMarkers ++;
        }
      
      // if the number of valid markers exceeds the 
      if (FilterByVesselness && nValidMarkers > nMarkers)
        {
        fExclude = true;
        }

      TransformType::Pointer transform = labelLineFilter->GetLineTransform();
      TransformType::MatrixType matrix = transform->GetMatrix();
      TransformType::OutputVectorType trans = transform->GetTranslation();

      point[0] = trans[0];
      point[1] = trans[1];
      point[2] = trans[2];
      norm[0]  = matrix[2][0];
      norm[1]  = matrix[2][1];
      norm[2]  = matrix[2][2];

      typedef LabelStatisticsType::LabelPixelType LabelPixelType;
      LabelPixelType labelValue = label;
      std::cout << "Detected line #"
                << label << ": "
                << "Exclude=" << fExclude
                << "; Point=("
                << point[0] << ", "
                << point[1] << ", "
                << point[2] << "); "
                << "Normal=("
                << norm[0] << ", "
                << norm[1] << ", "
                << norm[2] << "); "
                << "axisLength=("
                << axisLength[0] << ", "
                << axisLength[1] << ", "
                << axisLength[2] << "); "
                << "vessleness(min/max/mean/sigma)=("
                << labelStatistics->GetMinimum( labelValue ) << ", "
                << labelStatistics->GetMaximum( labelValue ) << ", "
                << labelStatistics->GetMean( labelValue ) << ", "
                << labelStatistics->GetSigma( labelValue ) << ")"
                << std::endl;
              
      if (fExclude)
        {
        changeMap[label] = 0;
        continue;
        }
      
      fixedPointSetContainer->InsertElement(pointId, point);
      pointId ++;
      fixedPointSetContainer->InsertElement(pointId, norm);
      pointId ++;

      Centroid[0] += point[0];
      Centroid[1] += point[1];
      Centroid[2] += point[2];
      
      }
    }
  
  double nPoints = (double)pointId/2.0;
  Centroid[0] /= nPoints;
  Centroid[1] /= nPoints;
  Centroid[2] /= nPoints;

  fixedPointSet->SetPoints( fixedPointSetContainer );

  // Generate a label map that only contains the identified marker
  // regions
  ChangeLabelFilter::Pointer changeLabel = ChangeLabelFilter::New();
  changeLabel->SetInput( RelabelFilter->GetOutput() );
  changeLabel->SetChangeMap( changeMap );
  writer->SetInput( changeLabel->GetOutput() );
  writer->SetUseCompression(1);

  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject &err)
    {
    std::cerr << err << std::endl;
    return EXIT_FAILURE ;
    }


  movingPointSet->SetPoints( movingPointSetContainer );

  typedef LineDistanceMetric::TransformType          TransformBaseType;
  typedef TransformBaseType::ParametersType         ParametersType;
  typedef TransformBaseType::JacobianType           JacobianType;

  LineDistanceMetric::Pointer metric = LineDistanceMetric::New();
  typedef itk::Euler3DTransform< double >      RegistrationTransformType;
  RegistrationTransformType::Pointer registrationTransform = RegistrationTransformType::New();

  // Optimizer Type
  typedef itk::LevenbergMarquardtOptimizer OptimizerType;
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  optimizer->SetUseCostFunctionGradient(false);


  // Registration Method
  typedef itk::PointSetToPointSetRegistrationMethod< PointSetType, PointSetType > RegistrationType;
  RegistrationType::Pointer   registration  = RegistrationType::New();
  //registration->SetNumberOfThreads(2);
  
  // Scale the translation components of the Transform in the Optimizer
  OptimizerType::ScalesType scales( registrationTransform->GetNumberOfParameters() );
  scales.Fill(0.01);
  //const double translationScale = 500;   // dynamic range of translations
  //const double rotationScale    = 5000;   // dynamic range of rotations
  //scales[0] = 1.0 / rotationScale;
  //scales[1] = 1.0 / rotationScale;
  //scales[2] = 1.0 / rotationScale;
  //scales[3] = 1.0 / translationScale; 
  //scales[4] = 1.0 / translationScale; 
  //scales[5] = 1.0 / translationScale;

  unsigned long   numberOfIterations =  1000;
  double          gradientTolerance  =  1e-5;   // convergence criterion
  double          valueTolerance     =  1e-5;   // convergence criterion
  double          epsilonFunction    =  1e-6;   // convergence criterion

  optimizer->SetScales( scales );
  optimizer->SetNumberOfIterations( numberOfIterations );
  optimizer->SetValueTolerance( valueTolerance );
  optimizer->SetGradientTolerance( gradientTolerance );
  optimizer->SetEpsilonFunction( epsilonFunction );


  // Start from an Identity transform (in a normal case, the user 
  // can probably provide a better guess than the identity...
  registrationTransform->SetIdentity();
  RegistrationTransformType::OutputVectorType vec;
  vec[0] = Centroid[0];
  vec[1] = Centroid[1];
  vec[2] = Centroid[2];
  registrationTransform->SetOffset(vec);

  registration->SetInitialTransformParameters( registrationTransform->GetParameters() );

  //------------------------------------------------------
  // Connect all the components required for Registration
  //------------------------------------------------------
  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetTransform(     registrationTransform     );
  registration->SetFixedPointSet( fixedPointSet );
  registration->SetMovingPointSet(   movingPointSet   );

  //Connect an observer
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  try 
    {
    registration->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    return EXIT_FAILURE;
    }

  // Convret Euler 3D transform to Rigid 3D Transform
  typedef itk::AffineTransform<double, 3> MarkerTransformType;
  //typedef typename itk::CenteredAffineTransform< double, 3 > MarkerTransformType;
  MarkerTransformType::Pointer transform = MarkerTransformType::New();

  transform->SetIdentity();
  transform->SetOffset(registrationTransform->GetOffset());
  transform->SetMatrix(registrationTransform->GetMatrix());
  
  std::cout << "Offset:   " << transform->GetOffset() << std::endl;
  std::cout << "Rotation: " << transform->GetMatrix() << std::endl;

  // Calculate Inverse Matrix to pass the transform to Slicer.
  MarkerTransformType::Pointer outputTransform = MarkerTransformType::New();
  transform->GetInverse(outputTransform);

  if (markerTransform != "")
    {
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer markerTransformWriter;
    markerTransformWriter= TransformWriterType::New();
    markerTransformWriter->SetFileName( markerTransform );
    markerTransformWriter->SetInput( outputTransform );
    try
      {
      markerTransformWriter->Update();
      }
    catch (itk::ExceptionObject &err)
      {
      std::cerr << err << std::endl;
      return EXIT_FAILURE ;
      }
    }

  return EXIT_SUCCESS;
}


} // end of anonymous namespace


int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType (inputVolume, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch (componentType)
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0));
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<char>(0));
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( argc, argv, static_cast<unsigned short>(0));
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0));
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( argc, argv, static_cast<unsigned int>(0));
        break;
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0));
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0));
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0));
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0));
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0));
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }

  catch( itk::ExceptionObject &excep)
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

/*=========================================================================

  Program:   3D Slicer Marker Registration CLI
  Module:    MarkerRegistration
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
#include "MarkerRegistrationCLP.h"

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

    std::cout << "Value = " << optimizer->GetCachedValue() << std::endl; 
    std::cout << "Position = "  << optimizer->GetCachedCurrentPosition();
    std::cout << std::endl << std::endl;

    }
   
};


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

  MovingPointSetContainer::Pointer movingPointContainer = MovingPointSetContainer::New();

  labelLineFilter->SetInput( RelabelFilter->GetOutput() );
  //typedef typename RelabelType::ObjectSizeInPhysicalUnitsContainerType PhysicalSizeContainerType; // for ITK v4
  typedef std::vector<float> PhysicalSizeContainerType;
  PhysicalSizeContainerType objectSize = RelabelFilter->GetSizeOfObjectsInPhysicalUnits();
  int nObjects = RelabelFilter->GetNumberOfObjects();

  // Setup a map for change label filter
  typedef itk::ChangeLabelImageFilter<OutputImageType, OutputImageType> ChangeLabelFilter;
  typedef typename ChangeLabelFilter::ChangeMapType ChangeMapType;
  //typename std::map< ChangeLabelFilter::InputPixelType, ChangeLabelFilter::OutputPixelType >  ChangeMapType;

  ChangeMapType changeMap;

  FixedPointSetType::Pointer fixedPointSet = FixedPointSetType::New();
  FixedPointSetContainer::Pointer fixedPointSetContainer = FixedPointSetContainer::New();
  unsigned int pointId = 0;
  for (int i= 0; i < nObjects; i ++) // Label 0 is background and skipped
    {
    // According to ITK's manual:
    // "Once all the objects are relabeled, the application can query the number of objects and
    //  the size of each object. Object sizes are returned in a vector. The size of the background
    //  is not calculated. So the size of object #1 is GetSizeOfObjectsInPixels()[0], the size of object
    //  #2 is GetSizeOfObjectsInPixels()[1], etc."
    float size = objectSize[i];
    int label = i + 1;
    if (size < minimumObjectSize || size > maximumObjectSize)
      {
      // Out of size criteria
      changeMap[label] = 0;
      }
    else
      {
      // TODO: It is probably a good idea to also check the length and thickness of
      // the segmented area using LabelToLineImageFilter::GetAxisLength().
      FixedPointType point;
      FixedPointType norm;
      labelLineFilter->SetLabel( label );
      labelLineFilter->Update();
      TransformType::Pointer transform = labelLineFilter->GetLineTransform();
      TransformType::MatrixType matrix = transform->GetMatrix();
      TransformType::OutputVectorType trans = transform->GetTranslation();
      point[0] = trans[0];
      point[1] = trans[1];
      point[2] = trans[2];
      norm[0]  = matrix[2][0];
      norm[1]  = matrix[2][1];
      norm[2]  = matrix[2][2];
      std::cerr << "Detected line #"
                << label
                << ": Point=("
                << point[0] << ", "
                << point[1] << ", "
                << point[2] << "); Normal=("
                << norm[0] << ", "
                << norm[1] << ", "
                << norm[2] << ")"
                << std::endl;
      fixedPointSetContainer->InsertElement(pointId, point);
      pointId ++;
      fixedPointSetContainer->InsertElement(pointId, norm);
      pointId ++;
      }
    }

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


  //
  // Create pointset of the Z-frame
  // TODO: this has to be configurable (e.g. XML)
  //
  MovingPointSetType::Pointer movingPointSet = MovingPointSetType::New();
  MovingPointSetContainer::Pointer movingPointSetContainer = MovingPointSetContainer::New();

  // int directions
  double R=-1.0; double L=1.0;
  double A=-1.0; double P=1.0;
  double S=1.0; double I=-1.0;
  pointId = 0;

  MovingPointType point0;
  MovingPointType norm0;
  point0[0] = L*30.0; point0[1] = P*30.0; point0[2] = S*0.0;
  norm0[0] = 0.0; norm0[1] = 0.0; norm0[2] = S*1.0;
  movingPointSetContainer->InsertElement(pointId, point0);
  pointId++;
  movingPointSetContainer->InsertElement(pointId, norm0);
  pointId++;

  MovingPointType point1;
  MovingPointType norm1;
  point1[0] = L*30.0; point1[1] = P*0.0; point1[2] = S*0.0;
  norm1[0] = 0.0; norm1[1] = A/sqrt(2.0); norm1[2] = S/sqrt(2.0);
  movingPointSetContainer->InsertElement(pointId, point1);
  pointId++;
  movingPointSetContainer->InsertElement(pointId, norm1);
  pointId++;

  MovingPointType point2;
  MovingPointType norm2;
  point2[0] = L*30.0; point2[1] = A*30.0; point2[2] = S*0.0;
  norm2[0] = 0.0; norm2[1] = 0.0; norm2[2] = S*1.0;
  movingPointSetContainer->InsertElement(pointId, point2);
  pointId++;
  movingPointSetContainer->InsertElement(pointId, norm2);
  pointId++;

  MovingPointType point3;
  MovingPointType norm3;
  point3[0] = L*0.0; point3[1] = A*30.0; point3[2] = S*0.0;
  norm3[0] = R/sqrt(2.0); norm3[1] = 0.0; norm3[2] = S/sqrt(2.0);
  movingPointSetContainer->InsertElement(pointId, point3);
  pointId++;
  movingPointSetContainer->InsertElement(pointId, norm3);
  pointId++;

  MovingPointType point4;
  MovingPointType norm4;
  point4[0] = R*30.0; point4[1] = A*30.0; point4[2] = S*0.0;
  norm4[0] = 0.0; norm4[1] = 0.0; norm4[2] = S*1.0;
  movingPointSetContainer->InsertElement(pointId, point4);
  pointId++;
  movingPointSetContainer->InsertElement(pointId, norm4);
  pointId++;

  MovingPointType point5;
  MovingPointType norm5;
  point5[0] = R*30.0; point5[1] = A*0.0; point5[2] = S*0.0;
  norm5[0] = 0.0; norm5[1] = P/2.0; norm5[2] = S/2.0;
  movingPointSetContainer->InsertElement(pointId, point5);
  pointId++;
  movingPointSetContainer->InsertElement(pointId, norm5);
  pointId++;

  MovingPointType point6;
  MovingPointType norm6;
  point6[0] = R*30.0; point6[1] = P*30.0; point6[2] = S*0.0;
  norm6[0] = 0.0; norm6[1] = 0.0; norm6[2] = S;
  movingPointSetContainer->InsertElement(pointId, point6);
  pointId++;
  movingPointSetContainer->InsertElement(pointId, norm6);
  pointId++;

  movingPointSet->SetPoints( movingPointSetContainer );

  //typedef itk::EuclideanDistanceLineMetric<PointSetType, PointSetType> LineDistanceMetric;
  typedef LineDistanceMetric::TransformType          TransformBaseType;
  typedef TransformBaseType::ParametersType         ParametersType;
  typedef TransformBaseType::JacobianType           JacobianType;

  LineDistanceMetric::Pointer metric = LineDistanceMetric::New();
  //typedef itk::Rigid3DTransform < double >   RegistrationTransformType;
  typedef itk::TranslationTransform< double, Dimension >      RegistrationTransformType;
  //typedef itk::Euler3DTransform< double >      RegistrationTransformType;
  RegistrationTransformType::Pointer registrationTransform = RegistrationTransformType::New();

  // Optimizer Type
  typedef itk::LevenbergMarquardtOptimizer OptimizerType;
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  optimizer->SetUseCostFunctionGradient(false);


  // Registration Method
  typedef itk::PointSetToPointSetRegistrationMethod< PointSetType, PointSetType > RegistrationType;
  RegistrationType::Pointer   registration  = RegistrationType::New();
  //registration->SetNumberOfThreads(1);
  
  // Scale the translation components of the Transform in the Optimizer
  OptimizerType::ScalesType scales( registrationTransform->GetNumberOfParameters() );
  scales.Fill( 0.1 );
  //scales.Fill( 1.0 );

  //unsigned long   numberOfIterations =  1000;
  //double          gradientTolerance  =  1e-5;    // convergence criterion
  //double          valueTolerance     =  1e-5;    // convergence criterion
  //double          epsilonFunction    =  1e-6;   // convergence criterion

  //const double translationScale = 10000.0;   // dynamic range of translations
  //const double rotationScale    =    1.0;   // dynamic range of rotations
  //
  //scales[0] = 1.0 / rotationScale;
  //scales[1] = 1.0 / rotationScale;
  //scales[2] = 1.0 / rotationScale;
  //scales[3] = 1.0 / translationScale; 
  //scales[4] = 1.0 / translationScale; 
  //scales[5] = 1.0 / translationScale;

  unsigned long   numberOfIterations =  5000;
  double          gradientTolerance  =  1e-4;   // convergence criterion
  double          valueTolerance     =  1e-4;   // convergence criterion
  double          epsilonFunction    =  1e-5;   // convergence criterion

  optimizer->SetScales( scales );
  optimizer->SetNumberOfIterations( numberOfIterations );
  optimizer->SetValueTolerance( valueTolerance );
  optimizer->SetGradientTolerance( gradientTolerance );
  optimizer->SetEpsilonFunction( epsilonFunction );


  // Start from an Identity transform (in a normal case, the user 
  // can probably provide a better guess than the identity...
  registrationTransform->SetIdentity();
  RegistrationTransformType::OutputVectorType vec;
  vec[0] = 0.0; vec[1] = 0.0; vec[2] = 0.0;
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
    registration->StartRegistration();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Solution = " << registrationTransform->GetParameters() << std::endl;

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

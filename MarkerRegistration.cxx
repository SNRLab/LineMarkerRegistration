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


#include "itkPluginUtilities.h"
#include "MarkerRegistrationCLP.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace {

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

  MovingPointSetType::Pointer movingPointSet = MovingPointSetType::New();
  MovingPointSetContainer::Pointer movingPointSetContainer = MovingPointSetContainer::New();
  //unsigned int pointId = 0;
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
      MovingPointType point;
      MovingPointType norm;
      labelLineFilter->SetLabel( label );
      labelLineFilter->Update();
      TransformType::Pointer transform = labelLineFilter->GetLineTransform();
      TransformType::MatrixType matrix = transform->GetMatrix();
      TransformType::OutputVectorType trans = transform->GetTranslation();
      point[0] = trans[0];
      point[1] = trans[1];
      point[2] = trans[2];
      norm[1]  = matrix[2][0];
      norm[2]  = matrix[2][0];
      norm[3]  = matrix[2][0];
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
      }
    }
  
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

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

  RelabelFilter->SetInput ( CCFilter->GetOutput() );
  RelabelFilter->SetMinimumObjectSize( minimumObjectSize );

  writer->SetInput( RelabelFilter->GetOutput() );
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

  // Detect lines from the label map
  labelLineFilter->SetInput( RelabelFilter->GetOutput() );

  OutputImageType::Pointer dummyImage = OutputImageType::New();

  // Line label 1
  typedef typename LabelLineFilterType::LineTransformType TransformType;
  labelLineFilter->SetLabel( 1 );
  TransformType::Pointer transform0 = labelLineFilter->GetLineTransform();
  labelLineFilter->Update();
  if (lineTransform0 != "")
    {
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer lineTransformWriter;
    lineTransformWriter= TransformWriterType::New();
    lineTransformWriter->SetFileName( lineTransform0 );
    lineTransformWriter->SetInput( transform0 );
    try
      {
      lineTransformWriter->Update();
      }
    catch (itk::ExceptionObject &err)
      {
      std::cerr << err << std::endl;
      return EXIT_FAILURE ;
      }
    }


  // Line label 2
  typedef typename LabelLineFilterType::LineTransformType TransformType;
  labelLineFilter->SetLabel( 2 );
  TransformType::Pointer transform1 = labelLineFilter->GetLineTransform();
  labelLineFilter->Update();
  if (lineTransform1 != "")
    {
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer lineTransformWriter;
    lineTransformWriter= TransformWriterType::New();
    lineTransformWriter->SetFileName( lineTransform1 );
    lineTransformWriter->SetInput( transform1 );
    try
      {
      lineTransformWriter->Update();
      }
    catch (itk::ExceptionObject &err)
      {
      std::cerr << err << std::endl;
      return EXIT_FAILURE ;
      }
    }

  // Line label 3
  typedef typename LabelLineFilterType::LineTransformType TransformType;
  labelLineFilter->SetLabel( 3 );
  TransformType::Pointer transform2 = labelLineFilter->GetLineTransform();
  labelLineFilter->Update();
  if (lineTransform2 != "")
    {
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer lineTransformWriter;
    lineTransformWriter= TransformWriterType::New();
    lineTransformWriter->SetFileName( lineTransform2 );
    lineTransformWriter->SetInput( transform2 );
    try
      {
      lineTransformWriter->Update();
      }
    catch (itk::ExceptionObject &err)
      {
      std::cerr << err << std::endl;
      return EXIT_FAILURE ;
      }
    }


  // Line label 4
  typedef typename LabelLineFilterType::LineTransformType TransformType;
  labelLineFilter->SetLabel( 4 );
  TransformType::Pointer transform3 = labelLineFilter->GetLineTransform();
  labelLineFilter->Update();
  if (lineTransform3 != "")
    {
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer lineTransformWriter;
    lineTransformWriter= TransformWriterType::New();
    lineTransformWriter->SetFileName( lineTransform3 );
    lineTransformWriter->SetInput( transform3 );
    try
      {
      lineTransformWriter->Update();
      }
    catch (itk::ExceptionObject &err)
      {
      std::cerr << err << std::endl;
      return EXIT_FAILURE ;
      }
    }

  // Line label 5
  typedef typename LabelLineFilterType::LineTransformType TransformType;
  labelLineFilter->SetLabel( 5 );
  TransformType::Pointer transform4 = labelLineFilter->GetLineTransform();
  labelLineFilter->Update();
  if (lineTransform4 != "")
    {
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer lineTransformWriter;
    lineTransformWriter= TransformWriterType::New();
    lineTransformWriter->SetFileName( lineTransform4 );
    lineTransformWriter->SetInput( transform4 );
    try
      {
      lineTransformWriter->Update();
      }
    catch (itk::ExceptionObject &err)
      {
      std::cerr << err << std::endl;
      return EXIT_FAILURE ;
      }
    }

  // Line label 6
  typedef typename LabelLineFilterType::LineTransformType TransformType;
  labelLineFilter->SetLabel( 6 );
  TransformType::Pointer transform5 = labelLineFilter->GetLineTransform();
  labelLineFilter->Update();
  if (lineTransform5 != "")
    {
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer lineTransformWriter;
    lineTransformWriter= TransformWriterType::New();
    lineTransformWriter->SetFileName( lineTransform5 );
    lineTransformWriter->SetInput( transform5 );
    try
      {
      lineTransformWriter->Update();
      }
    catch (itk::ExceptionObject &err)
      {
      std::cerr << err << std::endl;
      return EXIT_FAILURE ;
      }
    }

  // Line label 7
  typedef typename LabelLineFilterType::LineTransformType TransformType;
  labelLineFilter->SetLabel( 7 );
  TransformType::Pointer transform6 = labelLineFilter->GetLineTransform();
  labelLineFilter->Update();
  if (lineTransform6 != "")
    {
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer lineTransformWriter;
    lineTransformWriter= TransformWriterType::New();
    lineTransformWriter->SetFileName( lineTransform6 );
    lineTransformWriter->SetInput( transform6 );
    try
      {
      lineTransformWriter->Update();
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

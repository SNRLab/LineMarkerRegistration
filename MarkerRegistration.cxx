#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkLabelToNeedleImageFilter.h"

#include "itkOtsuThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkMinimumMaximumImageFilter.h"

//#include "itkMultiplyByConstantImageFilter.h"
#include "itkTransformFileWriter.h"

#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"


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
  //typedef   itk::ImageFileWriter< OutputImageType > WriterType;
  typedef   itk::ImageFileWriter< InternalImageType > WriterType;

  // Smoothing filter
  typedef   itk::SmoothingRecursiveGaussianImageFilter<
  InternalImageType, InternalImageType > SmoothingFilterType;
  
  // Line enhancement filter
  typedef   itk::Hessian3DToVesselnessMeasureImageFilter<
    InternalPixelType > LineFilterType;
  //typedef   itk::Hessian3DToNeedleImageFilter<
  //  InternalPixelType > LineFilterType;

  // Otsu Threshold Segmentation filter
  typedef   itk::OtsuThresholdImageFilter<
    InternalImageType, InternalImageType >  OtsuFilterType;
  typedef   itk::ConnectedComponentImageFilter<
    InternalImageType, OutputImageType >  CCFilterType;
  typedef   itk::RelabelComponentImageFilter<
    OutputImageType, OutputImageType > RelabelType;
  // Line detection filter
  typedef   itk::LabelToNeedleImageFilter<
  OutputImageType, OutputImageType > NeedleFilterType;

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

  typename ReaderType::Pointer reader = ReaderType::New();  
  typename WriterType::Pointer writer = WriterType::New();
  typename SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();

  /*
  typename HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
  */
  typename LineFilterType::Pointer lineFilter = LineFilterType::New();

  typename OtsuFilterType::Pointer OtsuFilter = OtsuFilterType::New();
  typename CCFilterType::Pointer CCFilter = CCFilterType::New();
  typename RelabelType::Pointer RelabelFilter = RelabelType::New();
  typename NeedleFilterType::Pointer needleFilter = NeedleFilterType::New();


  reader->SetFileName( inputVolume.c_str() );
  writer->SetFileName( outputVolume.c_str() );

  smoothing->SetInput( reader->GetOutput() );
  smoothing->SetSigma( static_cast< double >(sigma1) );

  lineFilter->SetAlpha1( static_cast< double >(alpha1));
  lineFilter->SetAlpha2( static_cast< double >(alpha2));
  
  MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter = MultiScaleEnhancementFilterType::New();
  multiScaleEnhancementFilter->SetInput(smoothing->GetOutput());
  multiScaleEnhancementFilter->SetSigmaMinimum(minsigma);
  multiScaleEnhancementFilter->SetSigmaMaximum(maxsigma);
  multiScaleEnhancementFilter->SetNumberOfSigmaSteps(stepsigma);
  multiScaleEnhancementFilter->SetHessianToMeasureFilter (lineFilter);

  //OtsuFilter->SetInput( multiScaleEnhancementFilter->GetOutput());
  //OtsuFilter->SetOutsideValue( 255 );
  //OtsuFilter->SetInsideValue(  0  );
  //OtsuFilter->SetNumberOfHistogramBins( numberOfBins );
  //CCFilter->SetInput (OtsuFilter->GetOutput());
  //CCFilter->FullyConnectedOff();
  //RelabelFilter->SetInput ( CCFilter->GetOutput() );
  //RelabelFilter->SetMinimumObjectSize( minimumObjectSize );
  //
  //needleFilter->SetInput( RelabelFilter->GetOutput() );
  //needleFilter->SetMinPrincipalAxisLength( static_cast< float >(minPrincipalAxisLength) );
  //needleFilter->SetAngleThreshold (static_cast< double >(anglethreshold) );
  //
  //// Set default orientation and closest point of the needle for detection
  //// Note that the parameter is passed in RAS coordinate system
  //// and must be converted to LPS coordinate system
  //needleFilter->SetNormal (static_cast< double >(-normal[0]),
  //                         static_cast< double >(-normal[1]),
  //                         static_cast< double >(normal[2]));
  //needleFilter->SetClosestPoint(static_cast< double >(-closestPoint[0]),
  //                              static_cast< double >(-closestPoint[1]),
  //                              static_cast< double >(closestPoint[2]));
  //

  //writer->SetInput( needleFilter->GetOutput() );
  writer->SetInput( multiScaleEnhancementFilter->GetOutput() );
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

  typedef typename NeedleFilterType::NeedleTransformType TransformType;
  TransformType::Pointer transform = needleFilter->GetNeedleTransform();

  if (needleTransform != "")
    {
    typedef itk::TransformFileWriter TransformWriterType;
    TransformWriterType::Pointer needleTransformWriter;
    needleTransformWriter= TransformWriterType::New();
    needleTransformWriter->SetFileName( needleTransform );
    needleTransformWriter->SetInput( transform );
    try
      {
      needleTransformWriter->Update();
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

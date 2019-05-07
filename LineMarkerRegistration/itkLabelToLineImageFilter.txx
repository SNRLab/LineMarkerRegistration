/*=========================================================================

  Program:   LineMarkerRegistration CLI for 3D Slicer
  Module:    itkLabelToLineImageFilter.h
  Language:  C++
  Contributor: Junichi Tokuda (BWH)

  This code is based on vtkImageToImageFilter.h in ITK.

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLabelToLineImageFilter_txx
#define __itkLabelToLineImageFilter_txx

#include "itkLabelToLineImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkListSample.h"
#include "itkCovarianceSampleFilter.h"
//#include "itkCovarianceCalculator.h"
#include "itkSymmetricEigenAnalysis.h"
#include "itkCrossHelper.h"

#include "vnl/vnl_math.h"

namespace itk
{

/**
 * Constructor
 */
template < typename  TInput, typename TOutput  >
LabelToLineImageFilter< TInput, TOutput >
::LabelToLineImageFilter()
{
  this->m_Label = 1;
  this->m_Normal[0] = 0.0;
  this->m_Normal[1] = 0.0;
  this->m_Normal[2] = 1.0;
}


template < typename  TInput, typename TOutput  >
void 
LabelToLineImageFilter< TInput, TOutput >
::GenerateData()
{
  itkDebugMacro(<< "LabelToLineImageFilter generating data ");
  
  typename InputImageType::ConstPointer input = this->GetInput();
  typename OutputImageType::Pointer output = this->GetOutput();

  ImageRegionConstIterator<InputImageType> it;
  it = ImageRegionConstIterator<InputImageType>( input, input->GetRequestedRegion() );
  ImageRegionIterator<OutputImageType> oit;
  this->AllocateOutputs();
  oit = ImageRegionIterator<OutputImageType>(output,
                                             output->GetRequestedRegion());

  output->FillBuffer(static_cast<OutputPixelType>(0));
  
  typedef itk::Statistics::ListSample< VectorType > PointListType;
  typedef std::map< InputPixelType, PointListType::Pointer > PointListMapType;

  //-- Create lists of points for selected label
  PointListType::Pointer plist = PointListType::New();

  oit.GoToBegin();
  it.GoToBegin();
  int found = 0;

  while (!it.IsAtEnd())
    {
    InputPixelType pix = it.Get();
    typename InputImageType::IndexType index = it.GetIndex();

    if (pix == this->m_Label)
      {
      found = 1;
      VectorType mv;
      typename InputImageType::PointType point;
      input->TransformIndexToPhysicalPoint (index, point);
      mv[0] = (double)point[0];
      mv[1] = (double)point[1];
      mv[2] = (double)point[2];
      plist->PushBack(mv);
      oit.Set(it.Get());
      }
    ++it;
    ++oit;
    }

  //-- For each label, perform principal component analysis
  VectorType lineNorm;      // Orientation normal vector of the line
  //VectorType lineTip;      // Tip of the line closest to the default point.
  VectorType lineCenterOfMass;      // Center of Mass

  InputPixelType pix = this->m_Label;
  PointListType::Pointer sample = plist; 
  //std::cout << "=== Label " << pix << "===" << std::endl;
  
  typedef itk::Statistics::CovarianceSampleFilter< PointListType > 
    CovarianceAlgorithmType;
  CovarianceAlgorithmType::Pointer covarianceAlgorithm = 
    CovarianceAlgorithmType::New();
  
  covarianceAlgorithm->SetInput( sample );
  covarianceAlgorithm->Update();
  
  //std::cout << "Sample covariance = " << std::endl ; 
  //std::cout << covarianceAlgorithm->GetCovarianceMatrix() << std::endl;
  
  CovarianceAlgorithmType::MeasurementVectorType meanVector;
  meanVector = covarianceAlgorithm->GetMean();
  //std::cout << "Sample mean = " << meanVector << std::endl ; 
  
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
  
  //std::cout << "EigenValues: " << eigenValues << std::endl;
  //std::cout << "EigenVectors (each row is an an eigen vector): " << std::endl;
  //std::cout << eigenMatrix << std::endl;

  // Check the direction of principal component
  VectorType principalVector = eigenMatrix[2];
  double ip = principalVector * m_Normal;
  
  // Calculate the line orientation vector.
  // If the default vector (m_Normal) and the principal
  // vector are oposit, flip the orientation.
  if (ip >= 0)
    {
    lineNorm = principalVector;
    }
  else
    {
    lineNorm = - principalVector;
    }
  
  lineNorm.Normalize();
  
  PointListType::Iterator iter = sample->Begin();
  VectorType vector;
  
  //// To detect the edge of the line artifact, calculate
  //// projections of the points in the line artifact, and
  //// find the farest from the center point (meanVector)

  VectorType axisVector[3];
  axisVector[0] = eigenMatrix[2]; // principalVector
  axisVector[1] = eigenMatrix[1];
  axisVector[2] = eigenMatrix[0];

  axisVector[0].Normalize();
  axisVector[1].Normalize();
  axisVector[2].Normalize();
    
  // Ranges of point distributions in 3 directions
  VectorType min;
  VectorType max;

  // Initialize range parameters
  vector = iter.GetMeasurementVector();
  for (int i = 0; i < 3; i ++)
    {
    min[i] = (vector-meanVector)*axisVector[i];
    max[i] = min[i];
    }
  
  while (iter != sample->End())
    {
    vector = iter.GetMeasurementVector();
    typename InputImageType::PointType point;
    typename OutputImageType::IndexType index;
    point[0] = vector[0];
    point[1] = vector[1];
    point[2] = vector[2];
    output->TransformPhysicalPointToIndex (point, index);
    output->SetPixel(index, pix);

    // Check the range in 3 directions
    for (int i = 0; i < 3; i ++)
      {
      double p = (vector-meanVector)*axisVector[i];
      if (p < min[i])
        {
        min[i] = p;
        }
      else if (p>max[i])
        {
        max[i] = p;
        }
      }
    
    ++ iter;
    }

  // Calculate the range in 3 directions
  m_AxisLength = max-min;
  
  //lineTip = meanVector + lineNorm * max;
  lineCenterOfMass = meanVector;
  
  // Output position and orientation of the line as Affine transform
  // (suppose an identitiy transform when the line orients (0, 0, 1) direction)
  m_LineTransform = LineTransformType::New();
  m_LineTransform->SetIdentity();

  if (found > 0)
    {
    VectorType nx;
    VectorType ny;
    nx[0] = 1.0; nx[1] = 0.0; nx[2] = 0.0;
    ny[0] = 0.0; ny[1] = 1.0; ny[2] = 0.0;
    typedef itk::CrossHelper< VectorType > CrossType;
    CrossType cross;
    VectorType t = cross(lineNorm, nx);
    VectorType s = cross(t, lineNorm);
    
    t.Normalize();
    s.Normalize();
    lineNorm.Normalize();
    
    LineTransformType::MatrixType matrix;
    matrix[0][0] = s[0];
    matrix[0][1] = s[1];
    matrix[0][2] = s[2];
    
    matrix[1][0] = t[0];
    matrix[1][1] = t[1];
    matrix[1][2] = t[2];
    
    matrix[2][0] = lineNorm[0];
    matrix[2][1] = lineNorm[1];
    matrix[2][2] = lineNorm[2];
    
    m_LineTransform->SetMatrix( matrix );
    m_LineTransform->Translate( lineCenterOfMass );
    }

}


template < typename  TInput, typename TOutput  >
void
LabelToLineImageFilter< TInput, TOutput >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}


} // end namespace itk
  
#endif

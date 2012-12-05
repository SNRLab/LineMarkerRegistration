/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkHessian3DToLineImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-04-25 12:27:26 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLabelToLineImageFilter_h
#define __itkLabelToLineImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkCenteredAffineTransform.h"

namespace itk
{
/** \class LabelToLineImageFilter
 * \brief 
 * \ingroup IntensityImageFilters TensorObjects
 *
 */
  
template < typename  TInput, typename TOutput  >
class ITK_EXPORT LabelToLineImageFilter : public
ImageToImageFilter< TInput, TOutput >
{
public:
  /** Standard class typedefs. */
  typedef LabelToLineImageFilter Self;
  typedef ImageToImageFilter<
          TInput,
          TOutput >                               Superclass;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  
  typedef typename Superclass::InputImageType            InputImageType;
  typedef typename Superclass::OutputImageType           OutputImageType;
  typedef typename InputImageType::PixelType             InputPixelType;
  typedef typename OutputImageType::PixelType            OutputPixelType;
  
  typedef typename itk::CenteredAffineTransform< float, 3 >      LineTransformType;

  typedef itk::Vector< double, 3 > VectorType;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( LabelToLineImageFilter, ImageToImageFilter );

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkSetMacro(Label, int);
  itkGetConstMacro(Label, int);

  // Returns axis length. axisLength[0] is the principal axis (along the line)
  inline void GetAxisLength(VectorType& axisLength)
  {
    axisLength[0] = m_AxisLength[0];
    axisLength[1] = m_AxisLength[1];
    axisLength[2] = m_AxisLength[2];
  }

  inline void SetNormal(double x, double y, double z)
  {
    m_Normal[0] = x; m_Normal[1] = y; m_Normal[2] = z;
  }

  inline void GetNormal(double& x, double& y, double& z)
  {
    x = m_Normal[0]; y = m_Normal[1]; z = m_Normal[2];
  }

  inline LineTransformType * GetLineTransform()
  {
    return m_LineTransform;
  }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(DoubleConvertibleToOutputCheck,
                  (Concept::Convertible<double, OutputPixelType>));
  /** End concept checking */
#endif


protected:
  LabelToLineImageFilter();
  ~LabelToLineImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Generate Data */
  void GenerateData( void );

private:
  LabelToLineImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  int m_Label;
  VectorType m_Normal;
  VectorType m_AxisLength;

  LineTransformType::Pointer m_LineTransform;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLabelToLineImageFilter.txx"
#endif
  
#endif

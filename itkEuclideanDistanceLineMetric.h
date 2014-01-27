/*=========================================================================

  Program:   LineMarkerRegistration CLI for 3D Slicer
  Module:    itkEuclideanDistanceLineMetric.h
  Language:  C++
  Contributor: Junichi Tokuda (BWH)

  This code is based on itkEuclideanDistanceMetric.h in ITK.

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkEuclideanDistanceLineMetric_h
#define __itkEuclideanDistanceLineMetric_h

#include "itkPointSetToPointSetMetric.h"
#include "itkCovariantVector.h"
#include "itkPoint.h"
#include "itkPointSet.h"
#include "itkImage.h"


// This is a quick and dirty implementation of a metric for line to line
// registration. Although this metric class bases PointSetPointSetMetric,
// it actually takes two sets of lines in the 3D space. The line sets
// are expressed by using PointSet type by storing point vectors and
// normal vectors as follows:
//
//  <point on line 1>, <normal vector for line 1>, <point on line 2>, <normal vector for line 2>, ... 
//
// For the future work, the following classes have to be replaced by LineSetToLineSet* classes
//    
//   itkPointSetToPointSetRegistrationMethod.h
//   itkPointSetToPointSetRegistrationMethod.txx
//   itkPointSetToPointSetMetric.h
//   itkPointSetToPointSetMetric.txx
//


namespace itk
{
/** \class EuclideanDistanceLineMetric
 * \brief Computes the minimum distance between a moving point-set
 *  and a fixed point-set. A vector of minimum closest point distance is
 *  created for each point in the moving point-set.
 *  No correspondance is needed.
 *  For speed consideration, the point-set with the minimum number of points
 *  should be used as the moving point-set.
 *  If the number of points is high, the possibility of setting a distance map
 *  should improve the speed of the closest point computation.
 *
 *  Reference: "A Method for Registration of 3-D Shapes",
 *             IEEE PAMI, Vol 14, No. 2, February 1992
 *
 * \ingroup RegistrationMetrics
 */
template < typename TFixedPointSet, typename TMovingPointSet, 
typename TDistanceMap = ::itk::Image<unsigned short, TMovingPointSet::PointDimension > >
class ITK_EXPORT EuclideanDistanceLineMetric : 
    public PointSetToPointSetMetric< TFixedPointSet, TMovingPointSet>
{
public:

  /** Standard class typedefs. */
  typedef EuclideanDistanceLineMetric                                Self;
  typedef PointSetToPointSetMetric<TFixedPointSet, TMovingPointSet >  Superclass;

  typedef SmartPointer<Self>         Pointer;
  typedef SmartPointer<const Self>   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
 
  /** Run-time type information (and related methods). */
  itkTypeMacro(EuclideanDistanceLineMetric, Object);
 
  /** Types transferred from the base class */
  typedef typename Superclass::TransformType              TransformType;
  typedef typename Superclass::TransformPointer           TransformPointer;
  typedef typename Superclass::TransformParametersType    TransformParametersType;
  typedef typename Superclass::TransformJacobianType      TransformJacobianType;

  typedef typename Superclass::MeasureType                MeasureType;
  typedef typename Superclass::DerivativeType             DerivativeType;
  typedef typename Superclass::FixedPointSetType          FixedPointSetType;
  typedef typename Superclass::MovingPointSetType         MovingPointSetType;
  typedef typename Superclass::FixedPointSetConstPointer  FixedPointSetConstPointer;
  typedef typename Superclass::MovingPointSetConstPointer MovingPointSetConstPointer;

  typedef typename Superclass::PointIterator              PointIterator;
  typedef typename Superclass::PointDataIterator          PointDataIterator;

  typedef TDistanceMap                                    DistanceMapType;
  typedef typename DistanceMapType::ConstPointer          DistanceMapPointer;


  /** Get the number of values */
  unsigned int GetNumberOfValues() const;

  double PointToLineDistanceSq( const typename Superclass::InputPointType& point,
                                const typename Superclass::InputPointType& lineBasePoint,
                                const typename Superclass::InputPointType& lineNormalVector) const;

  /** Get the derivatives of the match measure. */
  void GetDerivative( const TransformParametersType & parameters,
                      DerivativeType & Derivative ) const;

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue( const TransformParametersType & parameters ) const;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative( const TransformParametersType & parameters,
                              MeasureType& Value, DerivativeType& Derivative ) const;

  /** Set/Get the distance map */
  itkSetConstObjectMacro(DistanceMap,DistanceMapType);
  itkGetConstObjectMacro(DistanceMap,DistanceMapType);

  /** Set/Get if the distance should be squared. Default is true for computation speed */
  itkSetMacro(ComputeSquaredDistance,bool);
  itkGetConstMacro(ComputeSquaredDistance,bool);
  itkBooleanMacro(ComputeSquaredDistance);

protected:
  EuclideanDistanceLineMetric();
  virtual ~EuclideanDistanceLineMetric() {};

  /** PrintSelf funtion */
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  EuclideanDistanceLineMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  DistanceMapPointer m_DistanceMap;
  bool               m_ComputeSquaredDistance;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEuclideanDistanceLineMetric.txx"
#endif

#endif



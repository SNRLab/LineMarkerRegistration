/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkEuclideanDistancePointMetric.txx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkEuclideanDistanceLineMetric_txx
#define __itkEuclideanDistanceLineMetric_txx

#include "itkEuclideanDistanceLineMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

/** Constructor */
template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap> 
EuclideanDistanceLineMetric<TFixedPointSet,TMovingPointSet,TDistanceMap>
::EuclideanDistanceLineMetric()
{
  m_DistanceMap = 0;
  
  // when set to true it will be a bit faster, but it will result in minimizing
  // the sum of distances^4 instead of the sum of distances^2
  m_ComputeSquaredDistance = false; 
}

/** Return the number of values, i.e the number of points in the moving set */
template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
unsigned int
EuclideanDistanceLineMetric<TFixedPointSet,TMovingPointSet,TDistanceMap>  
::GetNumberOfValues() const
{
  MovingPointSetConstPointer movingPointSet = this->GetMovingPointSet();

  if( !movingPointSet ) 
    {
    itkExceptionMacro( << "Moving point set has not been assigned" );
    }

  return  movingPointSet->GetPoints()->Size();
}


/** Get the match Measure */
template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
typename EuclideanDistanceLineMetric<TFixedPointSet,TMovingPointSet,TDistanceMap>::MeasureType
EuclideanDistanceLineMetric<TFixedPointSet,TMovingPointSet,TDistanceMap>
::GetValue( const TransformParametersType & parameters ) const
{
  FixedPointSetConstPointer fixedPointSet = this->GetFixedPointSet();

  if( !fixedPointSet ) 
    {
    itkExceptionMacro( << "Fixed line set has not been assigned" );
    }
  else if (fixedPointSet->GetNumberOfPoints() % 2 != 0)
    {
    itkExceptionMacro( << "Fixed line set has wrong number of elements" );
    }

  
  MovingPointSetConstPointer movingPointSet = this->GetMovingPointSet();

  if( !movingPointSet ) 
    {
    itkExceptionMacro( << "Moving point set has not been assigned" );
    }
  else if (movingPointSet->GetNumberOfPoints() % 2 != 0)
    {
    itkExceptionMacro( << "Moving line set has wrong number of elements" );
    }


  PointIterator pointItr = movingPointSet->GetPoints()->Begin();
  PointIterator pointEnd = movingPointSet->GetPoints()->End();

  // length of measure is (number of point sets) / 2
  MeasureType measure;
  measure.set_size(movingPointSet->GetPoints()->Size()/2);

  this->SetTransformParameters( parameters );

  unsigned int identifier = 0;
  while( pointItr != pointEnd )
    {
    typename Superclass::InputPointType  inputPoint;
    inputPoint.CastFrom( pointItr.Value() );
    typename Superclass::OutputPointType transformedPoint = 
      this->m_Transform->TransformPoint( inputPoint );
    ++pointItr;

    typename Superclass::InputPointType  inputVector;
    inputVector.CastFrom( pointItr.Value() );
    typename Superclass::OutputPointType transformedVector = 
      this->m_Transform->TransformPoint( inputVector );
    ++pointItr;

    double minimunDistance = -1.0; // distance is initialized with value < 0
    bool closestPoint = false;
    if(!closestPoint)
      {
      // Find 2 points on the moving line sets
      typename Superclass::InputPointType::PixelType point0[PointDimensions];
      typename Superclass::InputPointType::PixelType point1[PointDimensions];
      for (int i = 0; i < PointDimensions)
        {
        // Assuming the length of the normal vector is 1
        point0[i] = transformedPoint->GetElement(i) - transformedVector->GetElement(i)/2.0;
        point1[i] = transformedPoint->GetElement(i) + transformedVector->GetElement(i)/2.0;
        }

      // Go trough the list of fixed point and find the closest distance
      PointIterator pointItr2 = fixedPointSet->GetPoints()->Begin();
      PointIterator pointEnd2 = fixedPointSet->GetPoints()->End();
    
      while( pointItr2 != pointEnd2 )
        {
        typename Superclass::InputPointType  lineBasePoint;
        lineBasePoint.CastFrom( pointItr2.Value() );
        ++pointItr2;
        
        typename Superclass::InputPointType  lineNormalVector;
        lineNormalVector.CastFrom( pointItr2.Value() );
        ++pointItr2;
        double sqdist0  = PointToLineDistance(point0, lineBasePoint, lineNormalVector);
        double sqdist1  = PointToLineDistance(point1, lineBasePoint, lineNormalVector);
        double rms_dist = vcl_sqrt((sqdist0+sqdist1)/2.0);
    
        if (minimumDistance < 0.0 || rms_dist < minimumDistance)
          {
          minimumDIstance = rms_dist;
          }
        }
      }
    measure.put(identifier,minimumDistance);
    
    ++pointItr;
    identifier++;
    }
  
  return measure;

}


template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
void
EuclideanDistanceLineMetric<TFixedPointSet,TMovingPointSet,TDistanceMap>
::PointToLineDistanceSq( typename InputPointType::PixelType *point,
                       typename InputPointType::PixelType *lineBasePoint,
                       typename InputPointType::PixelType *lineNormalVector) 
{
  // Calculate the distance between the point ('point') and the line represented
  // by the combination of the base point ('lineBasePoint') and the normal vector
  // ('lineNOrmalVector').

  double inner=0.0;
  for (int i = 0; i < PointDimensions)
    {
    inner += point[i]-lineBasePoint[i] * lineNormalVector[i];
    }
  
  double sqdistance=0.0;
  for (int i = 0; i < PointDimensions)
    {
    double elm = inner*lineNormalVector[i] + lineBasePoint[i] - point[i];
    sqdistance = elm*elm;
    }
  //double distance = vcl_sqrt(sqdistance);

  return sqdistance;
  
}


/** Get the Derivative Measure */
template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
void
EuclideanDistanceLineMetric<TFixedPointSet,TMovingPointSet,TDistanceMap>
::GetDerivative( const TransformParametersType & itkNotUsed(parameters),
                 DerivativeType & itkNotUsed(derivative) ) const
{

}

/** Get both the match Measure and theDerivative Measure  */
template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
void
EuclideanDistanceLineMetric<TFixedPointSet,TMovingPointSet,TDistanceMap>
::GetValueAndDerivative(const TransformParametersType & parameters, 
                        MeasureType & value, DerivativeType  & derivative) const
{
  value = this->GetValue(parameters);
  this->GetDerivative(parameters,derivative);
}

/** PrintSelf method */
template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>  
void
EuclideanDistanceLineMetric<TFixedPointSet,TMovingPointSet,TDistanceMap>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "DistanceMap: " << m_DistanceMap << std::endl;
  if(m_ComputeSquaredDistance)
    {
    os << indent << "m_ComputeSquaredDistance: True"<< std::endl;
    }
  else
    {
    os << indent << "m_ComputeSquaredDistance: False"<< std::endl;
    }
}

} // end namespace itk


#endif




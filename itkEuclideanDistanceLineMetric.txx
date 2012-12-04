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


    //double minimumDistance = NumericTraits<double>::max();
    //bool closestPoint = false;
    //
    //// Try to use the distance map to solve the closest point
    //if(m_DistanceMap)
    //  {
    //  // If the point is inside the distance map
    //  typename DistanceMapType::IndexType index;
    //  if(m_DistanceMap->TransformPhysicalPointToIndex(transformedPoint,index))
    //    {
    //    minimumDistance = m_DistanceMap->GetPixel(index);
    //    // In case the provided distance map was signed, 
    //    // we correct here the distance to take its absolute value.
    //    if( minimumDistance < 0.0 ) 
    //      {
    //      minimumDistance = -minimumDistance;
    //      }
    //    closestPoint = true;
    //    }
    //  }
    //
    //// if the closestPoint has not been found
    //if(!closestPoint)
    //  {
    //  // Go trough the list of fixed point and find the closest distance
    //  PointIterator pointItr2 = fixedPointSet->GetPoints()->Begin();
    //  PointIterator pointEnd2 = fixedPointSet->GetPoints()->End();
    //
    //  while( pointItr2 != pointEnd2 )
    //    {
    //    double dist = pointItr2.Value().SquaredEuclideanDistanceTo(transformedPoint);
    //
    //    if(!m_ComputeSquaredDistance)
    //      {
    //      dist = vcl_sqrt(dist);
    //      }
    //
    //    if(dist<minimumDistance)
    //      {
    //      minimumDistance = dist;
    //                }
    //    pointItr2++;
    //    }
    //  }

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
        double dist = pointItr2.Value().SquaredEuclideanDistanceTo(transformedPoint);
    
        if(!m_ComputeSquaredDistance)
          {
          dist = vcl_sqrt(dist);
          }
    
        if(dist<minimumDistance)
          {
          minimumDistance = dist;
                    }
        pointItr2++;
        }
      }
    measure.put(identifier,minimumDistance);

    ++pointItr;
    identifier++;
    }

  return measure;

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




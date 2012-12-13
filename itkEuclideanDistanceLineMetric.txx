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
  //measure.set_size(movingPointSet->GetPoints()->Size()/2);
  measure.set_size(movingPointSet->GetPoints()->Size());
  
  //if (movingPointSet->GetPoints()->Size() % 2 != 0)
  //  {
  //  return measure;
  //  }

  this->SetTransformParameters( parameters );
  unsigned int identifier = 0;

  std::vector<double> minimumDistanceVec;
  while( pointItr != pointEnd )
    {
    
    typename Superclass::InputPointType  inputPoint;

    // Base point
    inputPoint.CastFrom( pointItr.Value() );
    typename Superclass::OutputPointType transformedPoint = 
      this->m_Transform->TransformPoint( inputPoint );
    ++pointItr;

    // Normal vector
    inputPoint.CastFrom( pointItr.Value() );
    typename Superclass::OutputPointType transformedVector = 
      this->m_Transform->TransformPoint( inputPoint );
    ++pointItr;
    
    // Find 2 points on the moving line sets
    typename Superclass::InputPointType point0;
    typename Superclass::InputPointType point1;

    for (int i = 0; i < 3; i ++)
      {
      // Assuming the length of the normal vector is 1
      point0[i] = transformedPoint[i] - transformedVector[i]*30.0;
      point1[i] = transformedPoint[i] + transformedVector[i]*30.0;
      }
    // Go trough the list of fixed point and find the closest distance
    PointIterator pointItr2 = fixedPointSet->GetPoints()->Begin();
    PointIterator pointEnd2 = fixedPointSet->GetPoints()->End();
    
    double minimumDistance = NumericTraits<double>::max();
    double distance0 = 0.0;
    double distance1 = 0.0;

    while( pointItr2 != pointEnd2 )
      {
      typename Superclass::InputPointType  lineBasePoint;
      lineBasePoint.CastFrom( pointItr2.Value() );
      ++pointItr2;
      
      typename Superclass::InputPointType  lineNormalVector;
      lineNormalVector.CastFrom( pointItr2.Value() );
      ++pointItr2;

      //std::cerr << "lineBasePoint=" 
      //          << lineBasePoint[0] << ", "
      //          << lineBasePoint[1] << ", "
      //          << lineBasePoint[2] << std::endl;

      double sqdist0  = PointToLineDistanceSq(point0, lineBasePoint, lineNormalVector);
      double sqdist1  = PointToLineDistanceSq(point1, lineBasePoint, lineNormalVector);
      //double rms_dist = vcl_sqrt((sqdist0+sqdist1)/2.0);
      double rms_dist = vcl_sqrt(PointToLineDistanceSq(transformedPoint, lineBasePoint, lineNormalVector));
      double dist0  = vcl_sqrt(sqdist0);
      double dist1  = vcl_sqrt(sqdist1);
      if (rms_dist < minimumDistance)
        {
        minimumDistance = rms_dist;
        distance0 = dist0;
        distance1 = dist1;
        }
      }
    //minimumDistanceVec.push_back(minimumDistance);
    measure.put(identifier,minimumDistance);
    identifier++;
    measure.put(identifier,minimumDistance);
    identifier++;
    }
  
  //std::cerr << "minimumDistance = [";
  //std::vector<double>::iterator iter;
  //for (iter = minimumDistanceVec.begin(); iter != minimumDistanceVec.end(); iter ++)
  //  {
  //  std::cerr << *iter << ", ";
  //  }
  //std::cerr << "]" << std::endl;

  return measure;

}


template <class TFixedPointSet, class TMovingPointSet, class TDistanceMap>
double
EuclideanDistanceLineMetric<TFixedPointSet,TMovingPointSet,TDistanceMap>
::PointToLineDistanceSq( const typename Superclass::InputPointType& point,
                         const typename Superclass::InputPointType& lineBasePoint,
                         const typename Superclass::InputPointType& lineNormalVector) const
{
  // Calculate the distance between the point ('point') and the line represented
  // by the combination of the base point ('lineBasePoint') and the normal vector
  // ('lineNOrmalVector').

  double pointtopoint;
  double inner=0.0;
  for (int i = 0; i < 3; i ++)
    {
    inner += (point[i]-lineBasePoint[i]) * lineNormalVector[i];
    pointtopoint += (point[i]-lineBasePoint[i])*(point[i]-lineBasePoint[i]); //TEST
    }
  
  double sqdistance=0.0;
  for (int i = 0; i < 3; i ++)
    {
    double elm = inner*lineNormalVector[i] + lineBasePoint[i] - point[i];
    sqdistance += elm*elm;
    }
  //double distance = vcl_sqrt(sqdistance);
  //if (pointtopoint > 30*30)
  //  return pointtopoint;
  //else
  //  return sqdistance;
  return pointtopoint;
  
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




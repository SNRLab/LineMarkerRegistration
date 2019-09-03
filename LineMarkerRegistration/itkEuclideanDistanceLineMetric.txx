/*=========================================================================

  Program:   LineMarkerRegistration CLI for 3D Slicer
  Module:    itkEuclideanDistanceLineMetric.h
  Language:  C++
  Contributor: Junichi Tokuda (BWH)

  This code is based on itkEuclideanDistanceMetric.txx in ITK.

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

  m_LineMatchFlag = NULL;
  
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


  MovingPointIterator pointItr = movingPointSet->GetPoints()->Begin();
  MovingPointIterator pointEnd = movingPointSet->GetPoints()->End();


  MeasureType measure;
  int nMovingPoints = movingPointSet->GetPoints()->Size();
  measure.set_size( nMovingPoints );

  if (m_LineMatchFlag.IsNull())
    {
    return measure;
    }

  m_LineMatchFlag->Reserve( nMovingPoints );
  for (int i = 0; i < nMovingPoints; i ++)
    {
    m_LineMatchFlag->SetElement(i, true);
    }
  
  //if (movingPointSet->GetPoints()->Size() % 2 != 0)
  //  {
  //  return measure;
  //  }

  this->SetTransformParameters( parameters );
  unsigned int identifier = 0;

  // 
  typename TransformType::InputVectorType tmpInputVector;
  typename TransformType::OutputVectorType tmpOutputVector;
  
  std::map< int, int > closestLine;            // closestLine[fixedLineIndex] == movingLineIndex 
  std::map< int, double > closestLineDistance; // closestLineDistance[fixedLineIndex] == distance
  
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
    // NOTE: since InputPointType is used to present vector,
    //  we first convirt it to InputVectorType, transform, and
    //  convert to OutputPointType
    inputPoint.CastFrom( pointItr.Value() );
    tmpInputVector[0] = inputPoint[0];
    tmpInputVector[1] = inputPoint[1];
    tmpInputVector[2] = inputPoint[2];
    tmpOutputVector = this->m_Transform->TransformVector( tmpInputVector );

    typename Superclass::OutputPointType transformedVector;
    transformedVector[0] = tmpOutputVector[0];
    transformedVector[1] = tmpOutputVector[1];
    transformedVector[2] = tmpOutputVector[2];
    ++pointItr;
    
    // Find 2 points on the moving line sets
    typename Superclass::InputPointType point0;
    typename Superclass::InputPointType point1;

    for (int i = 0; i < 3; i ++)
      {
      // Assuming the length of the normal vector is 1
      point0[i] = transformedPoint[i] - transformedVector[i]*10.0;
      point1[i] = transformedPoint[i] + transformedVector[i]*10.0;
      }
    // Go trough the list of fixed point and find the closest distance
    FixedPointIterator pointItr2 = fixedPointSet->GetPoints()->Begin();
    FixedPointIterator pointEnd2 = fixedPointSet->GetPoints()->End();
    
    double minimumDistance = NumericTraits<double>::max();
    double distance0 = 0.0;
    double distance1 = 0.0;

    int closestLineIndex = -1;
    int lineIndex = 0;

    double penalty = 100.0;
    while( pointItr2 != pointEnd2 )
      {
      typename Superclass::InputPointType  lineBasePoint;
      lineBasePoint.CastFrom( pointItr2.Value() );
      ++pointItr2;
      
      typename Superclass::InputPointType  lineNormalVector;
      lineNormalVector.CastFrom( pointItr2.Value() );
      ++pointItr2;

      //double dist     = std::sqrt(PointToLineDistanceSq(transformedPoint, lineBasePoint, lineNormalVector));
      double sqdist0  = PointToLineDistanceSq(point0, lineBasePoint, lineNormalVector);
      double sqdist1  = PointToLineDistanceSq(point1, lineBasePoint, lineNormalVector);
      double dist0    = std::sqrt(sqdist0);
      double dist1    = std::sqrt(sqdist1);
      double dist     = std::sqrt((sqdist0+sqdist1)/2.0);
      if (dist < minimumDistance)
        {
        minimumDistance = dist;
        distance0 = dist0;
        distance1 = dist1;
        closestLineIndex = lineIndex;
        }
      lineIndex += 2;
      }
    
    if (closestLine.find(closestLineIndex) == closestLine.end())
      {
      closestLine[closestLineIndex] = identifier;
      closestLineDistance[closestLineIndex] = minimumDistance;
      }
    else if (closestLineDistance[closestLineIndex] > minimumDistance)
      {
      m_LineMatchFlag->SetElement(closestLine[closestLineIndex], false);
      m_LineMatchFlag->SetElement(closestLine[closestLineIndex]+1, false);

      closestLine[closestLineIndex] = identifier;
      closestLineDistance[closestLineIndex] = minimumDistance;
      }
    else
      {
      m_LineMatchFlag->SetElement(identifier, false);
      m_LineMatchFlag->SetElement(identifier+1, false);
      }
    measure.put(identifier,distance0);
    identifier++;
    measure.put(identifier,distance1);
    identifier++;
    }

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

  //double pointtopoint;
  double inner=0.0;
  for (int i = 0; i < 3; i ++)
    {
    inner += (point[i]-lineBasePoint[i]) * lineNormalVector[i];
    }
  double sqdistance=0.0;
  for (int i = 0; i < 3; i ++)
    {
    double elm = inner*lineNormalVector[i] + lineBasePoint[i] - point[i];
    sqdistance += elm*elm;
    }

  //pointtopoint = 
  //  (point[0]-lineBasePoint[0])*(point[0]-lineBasePoint[0])+ //TEST
  //  (point[1]-lineBasePoint[1])*(point[1]-lineBasePoint[1])+ //TEST
  //  (point[2]-lineBasePoint[2])*(point[2]-lineBasePoint[2]);
  
  //return pointtopoint;
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




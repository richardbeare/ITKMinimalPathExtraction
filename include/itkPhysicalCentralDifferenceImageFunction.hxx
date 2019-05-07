/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPhysicalCentralDifferenceImageFunction.hxx,v $
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkPhysicalCentralDifferenceImageFunction_hxx
#define itkPhysicalCentralDifferenceImageFunction_hxx

#include "itkPhysicalCentralDifferenceImageFunction.h"

namespace itk
{


/**
 * Constructor
 */
template <class TInputImage, class TCoordRep>
PhysicalCentralDifferenceImageFunction<TInputImage,TCoordRep>
::PhysicalCentralDifferenceImageFunction()
{
  m_Interpolator = InterpolateImageFunctionType::New();
  m_IllegalValue = NumericTraits< TCoordRep >::max();
}


/**
 *
 */
template <class TInputImage, class TCoordRep>
void
PhysicalCentralDifferenceImageFunction<TInputImage,TCoordRep>
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);
}


/**
 *
 */
template <class TInputImage, class TCoordRep>
typename PhysicalCentralDifferenceImageFunction<TInputImage,TCoordRep>::OutputType
PhysicalCentralDifferenceImageFunction<TInputImage,TCoordRep>
::Evaluate( const PointType& point ) const
{
  OutputType derivative;
  derivative.Fill( 0.0 );

  IndexType middle;
  m_Interpolator->GetInputImage()->TransformPhysicalPointToIndex(point, middle);

  for(unsigned int dim=0; dim<TInputImage::ImageDimension; dim++)
    {
    TCoordRep valueLeft, valueRight = m_IllegalValue;
    IndexType leftIndex(middle);
    IndexType rightIndex(middle);
    leftIndex[dim] -= 1;
    rightIndex[dim] += 1;

    // check the bounds
    // Get the left neighbor
    if (m_Interpolator->IsInsideBuffer(leftIndex))
      {
      const TCoordRep tmpLeft = m_Interpolator->GetInputImage()->GetPixel(leftIndex);
      if (tmpLeft < m_IllegalValue)
	{
	PointType pointLeft( point );
	pointLeft[dim] += -1 * Superclass::m_Image->GetSpacing()[dim];
	valueLeft = m_Interpolator->Evaluate( pointLeft );
	}
      }

    if (m_Interpolator->IsInsideBuffer(rightIndex))
      {
      const TCoordRep tmpRight = m_Interpolator->GetInputImage()->GetPixel(rightIndex);
      if (tmpRight < m_IllegalValue)
	{
	PointType pointRight( point );
	pointRight[dim] += 1 * Superclass::m_Image->GetSpacing()[dim];
	valueRight = m_Interpolator->Evaluate( pointRight );
	}
      }

    // Compute derivative
    if ((valueRight < m_IllegalValue) && (valueLeft < m_IllegalValue))
      {
      std::cout << "A " <<  middle << leftIndex << rightIndex << valueLeft << " " << valueRight << std::endl;
      derivative[dim] = (valueRight - valueLeft) *
	(0.5 / Superclass::m_Image->GetSpacing()[dim]);
      }
    else if ((valueRight > m_IllegalValue) && (valueLeft > m_IllegalValue))
      {
      std::cout << "B " <<  middle << leftIndex << rightIndex << valueLeft << " " << valueRight << std::endl;
      derivative[dim] = 0;
      }
    else if ((valueRight > m_IllegalValue) && (valueLeft > m_IllegalValue))
      {
      // This is where things get tricky. The right neighbour is
      // illegal, which means we need to use the left side only.
      // HOWEVER, if the right side is illegal, then our interpolated
      // estimate of our current location is probably influenced by
      // the illegal value. Thus we need to estimate the gradient from
      // the middle value and left
      TCoordRep middleValue =  m_Interpolator->GetInputImage()->GetPixel(middle);
      TCoordRep leftValueIdx = m_Interpolator->GetInputImage()->GetPixel(leftIndex);
      std::cout << "C " <<  middle << leftIndex << rightIndex << valueLeft << " " << valueRight << std::endl;
      derivative[dim] = (middleValue - leftValueIdx) /
	( Superclass::m_Image->GetSpacing()[dim]);
      }
    else if ((valueRight < m_IllegalValue) && (valueLeft > m_IllegalValue))
      {
      TCoordRep middleValue =  m_Interpolator->GetInputImage()->GetPixel(middle);
      TCoordRep rightValueIdx = m_Interpolator->GetInputImage()->GetPixel(rightIndex);
      std::cout << "D " <<  middle << leftIndex << rightIndex << valueLeft << " " << valueRight << std::endl;
      derivative[dim] = (rightValueIdx - middleValue) /
	( Superclass::m_Image->GetSpacing()[dim]);

      }
    }

  return ( derivative );
}


} // end namespace itk

#endif

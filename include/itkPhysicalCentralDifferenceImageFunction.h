/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPhysicalCentralDifferenceImageFunction.h,v $
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkPhysicalCentralDifferenceImageFunction_h
#define itkPhysicalCentralDifferenceImageFunction_h

#include "itkImageFunction.h"
#include "itkCovariantVector.h"
#include "itkImageBase.h"
#include "itkLinearInterpolateSelectedNeighborsImageFunction.h"

namespace itk
{

/**
 * \class PhysicalCentralDifferenceImageFunction
 * \brief Calculate the derivative by central differencing in physical space.
 *
 * This class is templated over the input image type and
 * the coordinate representation type (e.g. float or double).
 *
 * \author Dan Mueller, Queensland University of Technology, dan.muel[at]gmail.com
 *
 * \ingroup ImageFunctions
 *
 * \ingroup MinimalPathExtraction
 */

namespace Functor {  
template <typename TInput>
class ValidNeighbor
{
 public:
  ValidNeighbor() = default;
  ~ValidNeighbor() = default;
  bool operator!=(const ValidNeighbor &) const
  {
    return false;
  }

  bool operator==(const ValidNeighbor & other) const
  {
    return !( *this != other );
  }

  inline bool operator()(const TInput & A) const
  { return A < m_IllegalValue; }

 private:
  TInput m_IllegalValue = static_cast< TInput >( itk::NumericTraits< TInput >::max()/2 );
};
}
  
template <
  class TInputImage,
  class TCoordRep = float >
class ITK_EXPORT PhysicalCentralDifferenceImageFunction :
  public ImageFunction< TInputImage,
                        CovariantVector<double, TInputImage::ImageDimension>,
                        TCoordRep >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(PhysicalCentralDifferenceImageFunction);

  /** Dimension underlying input image. */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  /** Standard class type alias. */
  using Self = PhysicalCentralDifferenceImageFunction;
  using Superclass = ImageFunction<TInputImage,
                        CovariantVector<double,
                        itkGetStaticConstMacro(ImageDimension)>,
                        TCoordRep>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(PhysicalCentralDifferenceImageFunction, ImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** InputImageType type alias support. */
  using InputImageType = TInputImage;

  /** OutputType typdef support. */
  using OutputType = typename Superclass::OutputType;

  /** Index type alias support. */
  using IndexType = typename Superclass::IndexType;

  /** ContinuousIndex type alias support. */
  using ContinuousIndexType = typename Superclass::ContinuousIndexType;

  /** Point type alias support. */
  using PointType = typename Superclass::PointType;

  /** Linear interpolate function type alias. */
using InterpolateImageFunctionType = LinearInterpolateSelectedNeighborsImageFunction<TInputImage, TCoordRep, Functor::ValidNeighbor<typename TInputImage::PixelType> >;

  /** Set the input image.
   * \warning this method caches BufferedRegion information.
   * If the BufferedRegion has changed, user must call
   * SetInputImage again to update cached values. */
  void SetInputImage( const InputImageType * ptr ) override
    {
        this->Superclass::SetInputImage( ptr );
        if ( m_Interpolator.IsNotNull() )
          {
          m_Interpolator->SetInputImage( ptr );
          }
    }

  /** Evalulate the image derivative by central differencing at specified index.
   *
   *  No bounds checking is done. The point is assume to lie within the
   *  image buffer. ImageFunction::IsInsideBuffer() can be used to check
   *  bounds before calling this method. */
  OutputType EvaluateAtIndex( const IndexType& index ) const override
    {
    PointType point;
    m_Interpolator->GetInputImage()->TransformIndexToPhysicalPoint( index, point );
    return this->Evaluate( point );
    }

  /** Evalulate the image derivative by central differencing at specified
   *  continuous index.
   *
   *  No bounds checking is done. The point is assume to lie within the
   *  image buffer. ImageFunction::IsInsideBuffer() can be used to check
   *  bounds before calling this method. */
  OutputType EvaluateAtContinuousIndex(
    const ContinuousIndexType& cindex ) const override
    {
    PointType point;
    m_Interpolator->GetInputImage()->TransformContinuousIndexToPhysicalPoint( cindex, point );
    return this->Evaluate( point );
    }

  /** Evalulate the image derivative by central differencing at specified
   *  physical point.
   *
   *  No bounds checking is done. The point is assume to lie within the
   *  image buffer. ImageFunction::IsInsideBuffer() can be used to check
   *  bounds before calling this method. */
   OutputType Evaluate( const PointType& point ) const override;

protected:
  PhysicalCentralDifferenceImageFunction();
  ~PhysicalCentralDifferenceImageFunction() override{};
  void PrintSelf(std::ostream& os, Indent indent) const override;

private:
  typename InterpolateImageFunctionType::Pointer m_Interpolator;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkPhysicalCentralDifferenceImageFunction.hxx"
#endif

#endif

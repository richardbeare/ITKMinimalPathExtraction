/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSpeedFunctionPathInformation.h,v $
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkSpeedFunctionPathInformation_h
#define itkSpeedFunctionPathInformation_h

#include "itkLightObject.h"
#include "itkObjectFactory.h"

#include <vector>

namespace itk
{

/** \class SpeedFunctionPathInformation
 * \brief  PathInfo class for encapsulating information about a path
 * for a SpeedFunctionToPathFilter Object.
 * The points are stored as follows: end, start, way0, way1, ..., wayN.
 * Each element of a path can be a set of points - e.g. "end" may be
 * derived from a segmentation mask. Alternatively start and end may
 * be the left and right sides of a 2d image, producing a minimal
 * path across the image.
 * Fronts are propagated in reverse order: wayN, ..., way1, way0, start.
 * (NOTE: a front is never propagated from end).
 *
 * The user must provide at least one PathInfo object using
 * SpeedFunctionToPathFilter::AddPathInfo().
 * If multiple PathInfo objects are added,
 * multiple paths are extracted and saved to separate filter outputs.
 *
 * \author Dan Mueller,
 * Queensland University of Technology, dan.muel[at]gmail.com
 *
 * \sa LightObject
 *
 * \ingroup MinimalPathExtraction
 */
template <typename TPoint>
class SpeedFunctionPathInformation :
  public LightObject
{
public:
  /** Standard class typedefs. */
  typedef SpeedFunctionPathInformation Self;
  typedef LightObject                  Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  
  /** Run-time type information (and related methods). */
  itkTypeMacro(SpeedFunctionPathInformation, LightObject);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Some point typedefs. */
  typedef TPoint                   PointType;
  typedef std::vector< PointType > PointTypeVec;

  void ClearInfo();

  void SetStartPoint( const PointType & start );

  void SetEndPoint( const PointType & end );

  void AddWayPoint( const PointType & way );

  void SetCurrent( const PointType & current );

  /** Methods for adding extended path components **/
  void SetStartPoint( const PointTypeVec & start );

  void SetEndPoint( const PointTypeVec & end );

  void AddWayPoint( const PointTypeVec & way );

  void SetCurrent( const PointTypeVec & current );

  unsigned int GetNumberOfPoints( ) const;

  const PointTypeVec & GetStartPoint( ) const;

  const PointTypeVec & GetEndPoint( ) const;

  const PointTypeVec & GetWayPoint( SizeValueType i ) const;

  bool HasNextFront( ) const;

  const PointTypeVec & GetCurrentFrontAndAdvance( );

  const PointTypeVec & PeekCurrentFront( ) const;

  const PointTypeVec & PeekNextFront( ) const;

  const PointTypeVec & PeekPreviousFront( ) const;


protected:
  SpeedFunctionPathInformation( );
  virtual ~SpeedFunctionPathInformation( );
  virtual void PrintSelf( std::ostream& os, Indent indent ) const ITK_OVERRIDE;

  std::vector< PointTypeVec > m_Information;
  SizeValueType               m_Front;


private:
  ITK_DISALLOW_COPY_AND_ASSIGN(SpeedFunctionPathInformation);
};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSpeedFunctionPathInformation.hxx"
#endif

#endif

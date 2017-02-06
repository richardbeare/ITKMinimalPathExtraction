/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSpeedFunctionToPathFilter.hxx,v $
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef itkSpeedFunctionToPathFilter_hxx
#define itkSpeedFunctionToPathFilter_hxx

#include "vnl/vnl_math.h"
#include "itkSpeedFunctionToPathFilter.h"
#include "itkFastMarchingUpwindGradientImageFilter.h"


namespace itk
{

/**
 *
 */
template <class TInputImage, class TOutputPath>
SpeedFunctionToPathFilter<TInputImage,TOutputPath>
::SpeedFunctionToPathFilter()
{
  m_CurrentArrivalFunction = NULL;
}


/**
 *
 */
template <class TInputImage, class TOutputPath>
SpeedFunctionToPathFilter<TInputImage,TOutputPath>
::~SpeedFunctionToPathFilter()
{
}


/**
 *
 */
template<class TInputImage, class TOutputPath>
unsigned int
SpeedFunctionToPathFilter<TInputImage,TOutputPath>
::GetNumberOfPathsToExtract() const
{
  return m_Information.size();
}


/**
 *
 */
template<class TInputImage, class TOutputPath>
const typename SpeedFunctionToPathFilter<TInputImage,TOutputPath>::PointTypeVec &
SpeedFunctionToPathFilter<TInputImage,TOutputPath>
::GetNextEndPoint()
{
  return m_Information[Superclass::m_CurrentOutput]->GetEndPoint();
}


/**
 *
 */
template<class TInputImage, class TOutputPath>
typename SpeedFunctionToPathFilter<TInputImage,TOutputPath>::InputImageType *
SpeedFunctionToPathFilter<TInputImage,TOutputPath>
::ComputeArrivalFunction()
{
  // Get the speed image
  InputImagePointer speed =
    const_cast< InputImageType * >( this->GetInput() );

  // Set the fast marching method for computing the arrival function
  typedef FastMarchingUpwindGradientImageFilter< TInputImage, TInputImage >
    FastMarchingType;

  typedef typename FastMarchingType::NodeContainer NodeContainer;
  typedef typename FastMarchingType::NodeType      NodeType;
  typename FastMarchingType::Pointer marching = FastMarchingType::New();
  marching->SetInput( speed );
  marching->SetGenerateGradientImage( false );
  marching->SetTargetOffset( 2.0 * Superclass::m_TerminationValue );
  marching->SetTargetReachedModeToAllTargets( );

  // Add next and previous front sources as target points to
  // limit the front propagation to just the required zones
  PointTypeVec PrevFront = m_Information[Superclass::m_CurrentOutput]->PeekPreviousFront();
  PointTypeVec NextFront = m_Information[Superclass::m_CurrentOutput]->PeekNextFront();
  typedef std::vector < IndexType > IndexTypeVec;
  IndexTypeVec NextIndexVec(0);

  
  typename NodeContainer::Pointer targets = NodeContainer::New();
  targets->Initialize();

  for (typename PointTypeVec::iterator it = PrevFront.begin(); it != PrevFront.end(); it++)
    {
      IndexType indexTargetPrevious;
      NodeType nodeTargetPrevious;
      speed->TransformPhysicalPointToIndex( *it, indexTargetPrevious );
      nodeTargetPrevious.SetValue( 0.0 );
      nodeTargetPrevious.SetIndex( indexTargetPrevious );
      targets->InsertElement( 0, nodeTargetPrevious );
    }
  
  for (typename PointTypeVec::iterator it = NextFront.begin(); it != NextFront.end(); it++)
    {
      IndexType indexTargetNext;
      NodeType nodeTargetNext;
      speed->TransformPhysicalPointToIndex( *it, indexTargetNext );
      nodeTargetNext.SetValue( 0.0 );
      nodeTargetNext.SetIndex( indexTargetNext );
      targets->InsertElement( 1, nodeTargetNext );
      NextIndexVec.push_back(indexTargetNext);
    }
  marching->SetTargetPoints( targets );

  // Get the next Front source point and add as trial point
  typename NodeContainer::Pointer trial = NodeContainer::New();
  trial->Initialize();
  PointTypeVec CurrentFront =   m_Information[Superclass::m_CurrentOutput]->GetCurrentFrontAndAdvance();
  std::cout << PrevFront[0] << CurrentFront[0] << NextFront[0] << std::endl;
  std::cout << PrevFront.size() << " " << CurrentFront.size() << " " << NextFront.size() << std::endl;

  for (typename PointTypeVec::iterator it = CurrentFront.begin(); it != CurrentFront.end(); it++)
    {
      IndexType indexTrial;
      NodeType nodeTrial;
      speed->TransformPhysicalPointToIndex( *it, indexTrial );
      nodeTrial.SetValue( 0.0 );
      nodeTrial.SetIndex( indexTrial );
      trial->InsertElement( 0, nodeTrial );
    }
  marching->SetTrialPoints( trial );

  // Update the method and set the arrival function
  marching->UpdateLargestPossibleRegion( );
  m_CurrentArrivalFunction = marching->GetOutput( );
  m_CurrentArrivalFunction->DisconnectPipeline( );
  
  // Only the index with the minimum arrival time should stay in the "Next" point set
  // We've advanced, so we set the current to a single point
  if (NextFront.size() > 1) {
    InputImagePixelType MinTime = itk::NumericTraits<InputImagePixelType>::max();
    unsigned MinPos(0);
    for (unsigned idx = 0; idx < NextIndexVec.size(); ++idx)
      {
        InputImagePixelType V =  m_CurrentArrivalFunction->GetPixel(NextIndexVec[idx]);
        if (V < MinTime)
          {
            MinPos = idx;
            MinTime=V;
          }
      }
    std::cout << MinPos << " " << MinTime << " " << NextFront[MinPos] << std::endl;
    m_Information[Superclass::m_CurrentOutput]->SetCurrent(NextFront[MinPos]);
  }
  return m_CurrentArrivalFunction;
}


/**
 *
 */
template <class TInputImage, class TOutputPath>
void
SpeedFunctionToPathFilter<TInputImage,TOutputPath>
::GenerateData( void )
{
  // Get the speed function
  InputImagePointer speed =
    const_cast< InputImageType * >( this->GetInput() );
  if ( speed.IsNull() )
    {
    itkExceptionMacro( "Speed function image must be provided" );
    }

  // Ensure the user has added at least one path info object
  if ( m_Information.size() == 0 )
    {
    itkExceptionMacro( "No PathInfo objects: at least one must be added." );
    }

  // Extract the path
  Superclass::GenerateData( );
}


/**
 *
 */
template <class TInputImage, class TOutputPath>
void
SpeedFunctionToPathFilter<TInputImage,TOutputPath>
::Execute( const Object * object, const EventObject & event )
{
  (void) event;

  // Cast object to optmizer
  typename OptimizerType::Pointer optimizer = (OptimizerType*)
      dynamic_cast< const OptimizerType* >( object );
  if ( optimizer.IsNull() ) return;

  // Get current position and value
  typename OptimizerType::ParametersType currentParameters = optimizer->GetCurrentPosition();
  unsigned int lenParameters = currentParameters.GetSize();
  if ( lenParameters != InputImageDimension ) return;
  typename OptimizerType::MeasureType currentValue = optimizer->GetValue( currentParameters );

  // Convert parameters to point
  bool valid = false;
  unsigned int numparams = optimizer->GetCurrentPosition().GetSize();
  PointType point; point.Fill( 0.0 );
  for ( unsigned int i=0; i<numparams; i++ )
    {
    point[i] = optimizer->GetCurrentPosition()[i];
    valid = true;
    }
  if ( !valid ) return;

  // Check if we have reached the termination value
  if ( currentValue < Superclass::m_TerminationValue &&
       m_Information[Superclass::m_CurrentOutput]->HasNextFront() )
    {
    // We have terminated the current path segment,
    // but there are more fronts to propagate

    // TODO: The path has not actually reached the path point.
    //       Change the next front point to be the current point.

    // Update the arrival function and re-initialise the cost function
    Superclass::m_CostFunction->SetImage( this->ComputeArrivalFunction( ) );
    Superclass::m_CostFunction->Initialize( );
    }
  else if ( currentValue >= Superclass::m_TerminationValue )
    {
    // Convert point to continuous index
    InputImagePointer input = const_cast<InputImageType*>( this->GetInput() );
    ContinuousIndexType cindex;
    input->TransformPhysicalPointToContinuousIndex( point, cindex );

    // Add point as vertex in path
    OutputPathPointer output = this->GetOutput( Superclass::m_CurrentOutput );
    output->AddVertex( cindex );
    }
}


/**
 *
 */
template<class TInputImage, class TOutputPath>
void
SpeedFunctionToPathFilter<TInputImage,TOutputPath>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}


} // end namespace itk

#endif

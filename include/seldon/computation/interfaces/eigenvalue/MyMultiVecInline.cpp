#ifndef SELDON_FILE_MY_MULTIVEC_INLINE_CPP
#define SELDON_FILE_MY_MULTIVEC_INLINE_CPP

namespace Anasazi
{
  template<class ScalarType>
  inline int MyMultiVec<ScalarType>::GetVecLength () const
  {
    return(Length_);
  }
  
  
  template<class ScalarType>
  inline int MyMultiVec<ScalarType>::GetNumberVecs () const
  {
    return(NumberVecs_);
  }
  
  
  template<class ScalarType>
  inline ScalarType& MyMultiVec<ScalarType>
  ::operator()(const int i, const int j)
  {
    if (j < 0 || j >= NumberVecs_) throw(-1);
    if (i < 0 || i >= Length_) throw(-2);
    // 
    return(data_[j][i]);
  }
  
  
  template<class ScalarType>
  inline const ScalarType& MyMultiVec<ScalarType>
  ::operator()(const int i, const int j) const
  {
    if (j < 0 || j >= NumberVecs_) throw(-1);
    if (i < 0 || i >= Length_) throw(-2);
    // 
    return(data_[j][i]);
  }
  
  
  template<class ScalarType>
  inline ScalarType* MyMultiVec<ScalarType>::operator[](int v)
  {
    if (v < 0 || v >= NumberVecs_) throw(-1);
    return(data_[v]);
  }
  
  
  template<class ScalarType>
  inline ScalarType* MyMultiVec<ScalarType>::operator[](int v) const
  {
    return(data_[v]);
  }
  
}

#endif // MY_MULTIVECTOR_CPP

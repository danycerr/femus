#ifndef __dense_vectorM_h__
#define __dense_vectorM_h__

// #include "compare_types.h"
#include "dense_vector_baseM.h"

// C++ includes
#include <vector>
#include <cassert>
#include <cmath>



/**
 * Defines a dense vector for use in Finite Element-type computations.
 * This class is to basically compliment the \p DenseMatix class.  It
 * has additional capabilities over the \p std::vector that make it
 * useful for finite elements, particulary for systems of equations.
 *
 * @author Benjamin S. Kirk, 2003
 */ 

// ------------------------------------------------------------
// DenseVector class definition
class DenseVectorM : public DenseVectorBaseM
{
  private:

  /// The actual data values, stored as a 1D array.
  std::vector<double> _val;

public:

  /// Constructor.  Creates a dense vector of dimension \p n.
  explicit
  DenseVectorM(const unsigned int n=0);
  /// Copy-constructor.
  DenseVectorM (const DenseVectorM& other_vector);
  /// Copy-constructor, from a \p std::vector.
  DenseVectorM (const std::vector<double>& other_vector);
  /// Destructor.  Does nothing.     
  ~DenseVectorM() {}
  
  
   ///  * Resize the vector. Sets all elements to 0.
  void resize (const unsigned int n);
  /// Set every element in the vector to 0.
  virtual void zero();
  
  // ============================
  // RETURN FUNCTIONS
  // ===========================
  
   /// Puts the principal subvector of size \p sub_n
  void get_principal_subvector (unsigned int sub_n, DenseVectorM& dest) const;
  ///Access to the values array
  std::vector<double>& get_values() {return _val; }
  /// Access to the values array. 
  const std::vector<double>& get_values() const { return _val; }

  ///  * @returns the size of the vector.
  virtual unsigned int size() const { return _val.size(); }

  /// @returns the \p (i) element of the vector.
  double operator() (const unsigned int i) const;
  /// * @returns the \p (i,j) element of the vector as a writeable reference.
  double & operator() (const unsigned int i);

  ///  * @returns the \p (i) element of the vector.
  virtual double el(const unsigned int i) const { return (*this)(i); }
  ///  * @returns the \p (i) element of the vector as a writeable reference.
  virtual double & el(const unsigned int i)     { return (*this)(i); }
  
   // ============================
  // ALGEBRA
  // ===========================
  
  ///* Assignment operator.
  DenseVectorM& operator = (const DenseVectorM& other_vector);
  ///  * STL-like swap method
  void swap(DenseVectorM& other_vector);
  /// Multiplies every element in the vector by \p factor
  void scale (const double factor);
  
  /// * Multiplies every element in the vector by \p factor.
  DenseVectorM& operator*= (const double factor);
  
  ///  * Adds \p factor times \p vec to this vector. T += a * vec 
  void  add (const double a,  const DenseVectorM& vec);

  ///* Evaluate dot product with \p vect
  double dot (const DenseVectorM &vec) const;

  /// Tests if \p vec is exactly equal to this vector.
  bool operator== (const DenseVectorM &vec) const;
  /// Tests if \p vec is not exactly equal to this vector.
  bool operator!= (const DenseVectorM &vec) const;
  
  /// Adds \p vec to this vector.
  DenseVectorM& operator+= (const DenseVectorM &vec);
  
  /// Subtracts \p vec from this vector.
  DenseVectorM& operator-= (const DenseVectorM &vec);
  
  /// * @returns the minimum element in the vector.
  double min () const;

  /// @returns the maximum element in the vector.
  double max () const;

  /// @returns the \f$l_1\f$-norm of the vector, 
  double l1_norm () const;
  /// @returns the \f$l_2\f$-norm of the vector
  double l2_norm () const;
  /// @returns the maximum absolute value of the
  double linfty_norm () const;
};



// ------------------------------------------------------------
// DenseVector member functions

inline DenseVectorM::DenseVectorM(const unsigned int n) : _val (n, 0.){}


// ======================================================
inline DenseVectorM::DenseVectorM (const DenseVectorM& other_vector) :
  DenseVectorBaseM()
{
  const std::vector<double> &other_vals = other_vector.get_values();
  _val.clear();
  _val.reserve(other_vals.size());

  for (unsigned int i=0; i<other_vals.size(); i++)   _val.push_back(other_vals[i]);
}

// =================================
inline DenseVectorM::DenseVectorM (const std::vector<double>& other_vector) :
  _val(other_vector){  }


// =============================================
inline DenseVectorM& DenseVectorM::operator = (const DenseVectorM& other_vector)
{
  //  _val = other_vector._val;
  const std::vector<double> &other_vals = other_vector.get_values();
  _val.clear();  _val.reserve(other_vals.size());
  for (unsigned int i=0; i<other_vals.size(); i++)    _val.push_back(other_vals[i]);
  return *this;
}

// ======================================
inline void DenseVectorM::swap(DenseVectorM& other_vector){  _val.swap(other_vector._val);}
// =================================================
inline void DenseVectorM::resize(const unsigned int n){  _val.resize(n);  zero();}
// ==================
inline void DenseVectorM::zero(){  std::fill (_val.begin(),_val.end(),0.);}
// =================================
inline double DenseVectorM::operator () (const unsigned int i) const{assert(i < _val.size());  return _val[i]; }
// ====================================
inline double & DenseVectorM::operator () (const unsigned int i){   assert (i < _val.size()); return _val[i];}
    // ====================================  
inline void DenseVectorM::scale (const double factor){  for (unsigned int i=0; i<_val.size(); i++)    _val[i] *= factor;}
// =======================================
inline DenseVectorM& DenseVectorM::operator*= (const double factor){ this->scale(factor); return *this; }
// ======================================
inline void DenseVectorM::add (const double factor, const DenseVectorM& vec){
   assert(this->size() == vec.size());
  for (unsigned int i=0; i<this->size(); i++)    (*this)(i) += factor*vec(i);
}
// =========================================
inline double DenseVectorM::dot (const DenseVectorM& vec) const{ 
   assert(this->size() == vec.size());
  double val = 0.;
  for (unsigned int i=0; i<this->size(); i++)    val += (*this)(i)*vec(i);
  return val;
}
// ===========================================
inline bool DenseVectorM::operator== (const DenseVectorM& vec) const{
   assert (this->size() == vec.size());
  for (unsigned int i=0; i<this->size(); i++)    if ((*this)(i) != vec(i))
      return false;

  return true;
}

// ===============================================
inline bool DenseVectorM::operator!= (const DenseVectorM& vec) const
{
   assert (this->size() == vec.size());
  for (unsigned int i=0; i<this->size(); i++)    if ((*this)(i) != vec(i))
      return true;
  return false;
}

// ==================================
inline DenseVectorM& DenseVectorM::operator+= (const DenseVectorM& vec){
   assert (this->size() == vec.size());
  for (unsigned int i=0; i<this->size(); i++)    (*this)(i) += vec(i);
  return *this;
}

// ==============================
inline DenseVectorM& DenseVectorM::operator-= (const DenseVectorM& vec)
{
   assert (this->size() == vec.size());
  for (unsigned int i=0; i<this->size(); i++)    (*this)(i) -= vec(i);
  return *this;
}

// ===================================
inline double DenseVectorM::min () const
{
   assert (this->size());
  double my_min = (*this)(0);
  for (unsigned int i=1; i!=this->size(); i++)    {
      double current = (*this)(i);
      my_min = (my_min < current? my_min : current);
    }
  return my_min;
}

// ==================================
inline double DenseVectorM::max () const
{
   assert (this->size());  double my_max = (*this)(0);
  for (unsigned int i=1; i!=this->size(); i++){
      double current = (*this)(i);
      my_max = (my_max > current? my_max : current);
    }
  return my_max;
}

// =======================
inline double DenseVectorM::l1_norm () const
{
  double my_norm = 0.;
  for (unsigned int i=0; i!=this->size(); i++)  my_norm += fabs((*this)(i));
  return my_norm;
}

// ==================================================
inline double DenseVectorM::l2_norm () const
{
  double my_norm = 0.;
  for (unsigned int i=0; i!=this->size(); i++)   my_norm +=((*this)(i))*((*this)(i));
  return sqrt(my_norm);
}

// ======================================================
inline double DenseVectorM::linfty_norm () const{
  if (!this->size())    return 0.;
  double my_norm = fabs((*this)(0));

  for (unsigned int i=1; i!=this->size(); i++)    {
      double current = fabs((*this)(i));
      my_norm = (my_norm > current? my_norm : current);
    }
  return sqrt(my_norm);
}

// ==================================================================
inline void DenseVectorM::get_principal_subvector (unsigned int sub_n,
                                              DenseVectorM& dest) const
{
   assert( sub_n <= this->size() );
  dest.resize(sub_n);
  for(unsigned int i=0; i<sub_n; i++)  dest(i) = (*this)(i);
}

#endif // #ifndef __dense_vector_h__


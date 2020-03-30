#ifndef SELDON_FILE_VIRTUAL_EIGENVALUE_SOLVER_HXX

namespace Seldon
{
  
  //! Base class to solve an eigenvalue problem
  /*!
    Resolution of a standard eigenvalue problem : K x = lambda x
    or a generalized eigenvalue problem K x = lambda M x
    M is called mass matrix, K stiffness matrix, lambda is the eigenvalue
    and x the eigenvector.
    
    This class should not be instantiated directly, but rather derived classes
    like DenseEigenProblem, SparseEigenProblem, MatrixFreeEigenProblem.
  */
  template<class T>
  class EigenProblem_Base
  {
  public :
    //! several available modes to find eigenvalues (Arpack)
    /*!
      REGULAR_MODE : Regular mode
                 standard problem => no linear system to solve
                 generalized problem => M^-1 K x = lambda x (inverse of M required)
      SHIFTED_MODE : Shifted mode 
                 standard problem => (K - sigma I)^-1 X = lambda X
                 generalized problem => (K - sigma M)^-1 M X = lambda X
      BUCKLING_MODE : Buckling mode  (real symmetric problem)
                 generalized problem => (K - sigma M)^-1 K X = lambda X
      CAYLEY_MODE : Cayley mode (real symmetric problem)
                 generalized problem => (K - sigma M)^-1 (K + sigma M) X = lambda X
      INVERT_MODE : Shifted mode on matrix M^-1 K      
      IMAG_SHIFTED_MODE : using Imag( (K - sigma M)^-1 M)
                     instead of Real( (K - sigma M)^-1 M)
                         mode 4 in Arpack (dnaupd.f)
    */
    enum {REGULAR_MODE, SHIFTED_MODE, IMAG_SHIFTED_MODE, INVERT_MODE,
          BUCKLING_MODE, CAYLEY_MODE};
    
    //! parts of the spectrum (near from 0, at infinity or around a given value)
    /*!
      SMALL_EIGENVALUES : seeking eigenvalues near 0
      LARGE_EIGENVALUES : seeking largest eigenvalues
      CENTERED_EIGENVALUES : seeking eigenvalues near the shift sigma
    */
    enum {SMALL_EIGENVALUES, LARGE_EIGENVALUES, CENTERED_EIGENVALUES};

    //! different sorting strategies
    enum {SORTED_REAL, SORTED_IMAG, SORTED_MODULUS};

    //! real part, imaginary part or complex solution
    enum { REAL_PART, IMAG_PART, COMPLEX_PART};
    
    //! different solvers (Anasazi)
    /*!
      SOLVER_LOBPCG : Locally Optimal Block Preconditioned Conjugate Gradient
      SOLVER_BKS : Block Krylov Schur
      SOLVER_BD : Block Davidson
    */
    enum {SOLVER_LOBPCG, SOLVER_BKS, SOLVER_BD};
    
    //! orthogonalization managers (Anasazi)
    enum {ORTHO_DGKS, ORTHO_SVQB};
        
  protected :
    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;
    
    //! mode used to find eigenvalues (regular, shifted, Cayley, etc)
    int eigenvalue_computation_mode;
    
    //! number of eigenvalues to be computed
    int nb_eigenvalues_wanted; 
    
    //! additional number of eigenvalues
    /*! Sometimes Arpack finds more converged eigenvalues than asked
      it is needed to store these eigenvalues and eigenvalues
      to avoid segmentation fault */
    int nb_add_eigenvalues;

    //! which spectrum ? Near from Zero ? Near from Infinity ? or near from a value ?
    int type_spectrum_wanted;
    
    //! large eigenvalues because of their real part, imaginary part or magnitude ?
    int type_sort_eigenvalues;
    
    //! if true, the generalized problem is reduced to a standard problem
    /*!
      If matrix M is symmetric definite positive, one may compute Cholesky
      factorisation of M = L L^T, and find eigenvalues of the standard problem :
      L^-1 K L^-T x = lambda x
     */    
    bool use_cholesky;
    
    //! if true, the generalized problem is reduced to a standard problem
    /*!
      if M is diagonal, one can seek eigenvalues of the 
      standard problem M^-1/2  K  M^-1/2  x = lambda x
    */
    bool diagonal_mass;
      
    //! threshold for Arpack's iterative process
    Treal stopping_criterion;
    
    //! Maximal number of iterations
    int nb_maximum_iterations;
    
    //! number of matrix-vector products
    int nb_prod;
    
    //! size of the problem
    int n_;
    
    //! shift sigma (if type_spectrum = centered_eigenvalues)
    T shift, shift_imag;
    
    //! number of Arnoldi vectors
    int nb_arnoldi_vectors;
    
    //! if true nb_arnoldi_vectors is automatically computed
    bool automatic_selection_arnoldi_vectors;
    
    int print_level; //!< print level
    
    //! if true consider Real( (a M + bK)^-1) or Imag( (a M + b K)^-1 )
    //! or the whole system, a and/or b being complex
    bool complex_system; int selected_part;    
    
#ifdef SELDON_WITH_MPI
    //! communicator used to compute eigenvalues
    MPI::Comm* comm;
#endif

    //! which solver ?
    int type_solver;
    
    
    /************************** 
     * Parameters for Anasazi *
     **************************/
    
    //! orthogonalization manager
    int ortho_manager;
    
    //! number of blocks for blocked solvers
    int nb_blocks;
    
    //! restart parameter for blocked solvers
    int restart_number;
    
    /************************
     * Parameters for Feast *
     ************************/
    
    //! interval where eigenvalues are searched
    Treal emin_interval, emax_interval;
    
  public :

    EigenProblem_Base();
    virtual ~EigenProblem_Base();
    
    // initialization
    void Init(int n);
    
    // basic functions
    int GetComputationalMode() const;
    void SetComputationalMode(int mode);
    
#ifdef SELDON_WITH_MPI
    void SetCommunicator(MPI::Comm& comm_);
    MPI::Comm& GetCommunicator();
#endif
    
    int GetNbAskedEigenvalues() const;
    void SetNbAskedEigenvalues(int n);

    int GetNbAdditionalEigenvalues() const;
    void SetNbAdditionalEigenvalues(int n);
    
    int GetNbBlocks() const;
    void SetNbBlocks(int);
    int GetNbMaximumRestarts() const;
    void SetNbMaximumRestarts(int);
    int GetOrthoManager() const;

    int GetEigensolverType() const;
    void SetEigensolverType(int type);
    
    int GetTypeSpectrum() const;
    int GetTypeSorting() const;
    T GetShiftValue() const;
    T GetImagShiftValue() const;
    void SetShiftValue(const T&);
    void SetImagShiftValue(const T&);

    void GetComplexShift(const Treal&, const Treal&, Tcplx&) const;
    void GetComplexShift(const Tcplx&, const Tcplx&, Tcplx&) const;

    void SetTypeSpectrum(int type, const T& val,
                         int type_sort = SORTED_MODULUS);
    
    void SetTypeSpectrum(int type, const complex<T>& val,
                         int type_sort = SORTED_MODULUS);

    Treal GetLowerBoundInterval() const;
    Treal GetUpperBoundInterval() const;

    void SetIntervalSpectrum(Treal, Treal);
            
    void SetCholeskyFactoForMass(bool chol = true);
    bool UseCholeskyFactoForMass() const;
    
    void SetDiagonalMass(bool diag = true);
    bool DiagonalMass() const;
    
    void SetStoppingCriterion(Treal eps);
    Treal GetStoppingCriterion() const;
    void SetNbMaximumIterations(int n);
    int GetNbMaximumIterations() const;
    
    int GetNbMatrixVectorProducts() const;
    
    int GetNbArnoldiVectors() const;
    void SetNbArnoldiVectors(int n);
    
    int GetM() const;
    int GetN() const;
    
    int GetPrintLevel() const;
    void SetPrintLevel(int lvl);
    
    void IncrementProdMatVect();
    
    void PrintErrorInit() const;
    virtual bool IsSymmetricProblem() const = 0;
    virtual bool IsHermitianProblem() const = 0;
    
    // mass matrix stuff
    virtual void ComputeDiagonalMass() = 0;
    virtual void FactorizeDiagonalMass() = 0;
    virtual void GetSqrtDiagonal(Vector<T>&) = 0;
    
    virtual void MltInvSqrtDiagonalMass(Vector<Treal>& X) = 0;
    virtual void MltSqrtDiagonalMass(Vector<Treal>& X) = 0;

    virtual void MltInvSqrtDiagonalMass(Vector<Tcplx>& X) = 0;
    virtual void MltSqrtDiagonalMass(Vector<Tcplx>& X) = 0;
    
    virtual void ComputeMassForCholesky();
    
    virtual void ComputeMassMatrix();
    virtual void MltMass(const Vector<T>& X, Vector<T>& Y) = 0;

    // stiffness matrix stuff
    virtual void ComputeStiffnessMatrix();    
    virtual void ComputeStiffnessMatrix(const T& a, const T& b);
    
    virtual void MltStiffness(const Vector<T>& X, Vector<T>& Y) = 0;
    virtual void MltStiffness(const T& a, const T& b,
			      const Vector<T>& X, Vector<T>& Y) = 0;
    
    // functions to overload (factorisation of mass and/or stiffness matrix)
    virtual void ComputeAndFactorizeStiffnessMatrix(const Treal& a, const Treal& b,
						    int real_p = COMPLEX_PART);
    
    virtual void ComputeAndFactorizeStiffnessMatrix(const Tcplx& a, const Tcplx& b,
						    int real_p = COMPLEX_PART);
    
    virtual void ComputeSolution(const Vector<Treal>& X, Vector<Treal>& Y);
    virtual void ComputeSolution(const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    virtual void ComputeSolution(const SeldonTranspose&,
				 const Vector<Treal>& X, Vector<Treal>& Y);

    virtual void ComputeSolution(const SeldonTranspose&,
				 const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    virtual void FactorizeCholeskyMass();
    
    virtual void MltCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X);
    virtual void MltCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X);
    
    virtual void SolveCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X);
    virtual void SolveCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X);
    
    virtual void Clear();
    
    // modification of eigenvectors to take into account 
    // the use of matrix D^-1/2 K D^-1/2 or L^-1 K L^-T instead of K
    // => eigenvectors are recovered by multiplying them by matrix D^1/2 or by L^T
    template<class T0, class Prop, class Storage>
    friend void ApplyScalingEigenvec(EigenProblem_Base<T0>& var,
                                     Vector<T0>& eigen_values,
                                     Vector<T0>& lambda_imag,
                                     Matrix<T0, Prop, Storage>& eigen_vectors,
                                     const T0& shiftr, const T0& shifti);

    template<class T0, class Prop, class Storage>
    friend void ApplyScalingEigenvec(EigenProblem_Base<complex<T0> >& var,
                                     Vector<complex<T0> >& eigen_values,
                                     Vector<complex<T0> >& lambda_imag,
                                     Matrix<complex<T0>, Prop, Storage>& eigen_vectors,
                                     const complex<T0>& shiftr,
                                     const complex<T0>& shifti);
    
  };    

  
  //! base class for eigenvalue problems
  template<class T, class StiffValue = T,
	   class MassValue = typename ClassComplexType<T>::Treal>  
  class VirtualEigenProblem : public EigenProblem_Base<T>
  {
  protected:
    //! mass matrix
    VirtualMatrix<MassValue>* Mh;

    //! stiffness matrix
    VirtualMatrix<StiffValue>* Kh;

    //! diagonal D^1/2 if the mass matrix is diagonal positive
    Vector<MassValue> sqrt_diagonal_mass;

    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;

  public:    
    VirtualEigenProblem();
    
    void InitMatrix(VirtualMatrix<StiffValue>& K);
    void InitMatrix(VirtualMatrix<StiffValue>& K, VirtualMatrix<MassValue>& M);

    void SetMatrix(VirtualMatrix<StiffValue>& K, VirtualMatrix<MassValue>& M);
    
    void SetTypeSpectrum(int type, const T& val,
                         int type_sort = EigenProblem_Base<T>::SORTED_MODULUS);

    void SetTypeSpectrum(int type, const complex<T>& val,
                         int type_sort = EigenProblem_Base<T>::SORTED_MODULUS);

    bool IsSymmetricProblem() const;
    bool IsHermitianProblem() const;

    void ComputeDiagonalMass();
    void FactorizeDiagonalMass();
    void GetSqrtDiagonal(Vector<T>&);
    
    void MltInvSqrtDiagonalMass(Vector<Treal>& X);
    void MltSqrtDiagonalMass(Vector<Treal>& X);

    void MltInvSqrtDiagonalMass(Vector<Tcplx>& X);
    void MltSqrtDiagonalMass(Vector<Tcplx>& X);

    void MltMass(const Vector<T>& X, Vector<T>& Y);    
    
    void MltStiffness(const Vector<T>& X, Vector<T>& Y);
    void MltStiffness(const T& coef_mass, const T& coef_stiff,
		      const Vector<T>& X, Vector<T>& Y);
    void Clear();
    
  };
  
  
  // sorting eigenvalues
  template<class T, class Storage1, class Storage2>
  void SortEigenvalues(Vector<T>& lambda_r, Vector<T>& lambda_i,
                       Matrix<T, General, Storage1>& eigen_old,
                       Matrix<T, General, Storage2>& eigen_new,
                       int type_spectrum, int type_sort,
                       const T& shift_r, const T& shift_i);
  
  // sorting eigenvalues
  template<class T, class Storage1, class Storage2>
  void SortEigenvalues(Vector<complex<T> >& lambda_r,
                       Vector<complex<T> >& lambda_i,
                       Matrix<complex<T>, General, Storage1>& eigen_old,
                       Matrix<complex<T>, General, Storage2>& eigen_new,
                       int type_spectrum, int type_sort,
                       const complex<T>& shift_r, const complex<T>& shift_i);
  
  
  //! computation of a few eigenvalues for dense matrices
  template<class T, class Tstiff, class Prop, class Storage,
           class Tmass = typename ClassComplexType<T>::Treal,
	   class PropM = Symmetric, class StorageM = RowSymPacked>
  class DenseEigenProblem : public VirtualEigenProblem<T, Tstiff, Tmass>
  {
  protected :
    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;

    //! LU factorisation of a real matrix
    Matrix<Treal, Prop, Storage> mat_lu_real;

    //! LU factorisation of a complex matrix
    Matrix<Tcplx, Prop, Storage> mat_lu_cplx;

    //! mass matrix
    Matrix<Tmass, PropM, StorageM>* Mh;

    //! stiffness matrix
    Matrix<Tstiff, Prop, Storage>* Kh;
    
    //! pivot used by the LU factorisation
    Vector<int> pivot;
    
    //! Cholesky factorisation of mass matrix
    Matrix<Treal, PropM, StorageM> mat_chol;
    
    //! temporary vectors for Cholesky
    Vector<Treal> Xchol_real, Xchol_imag;
    
  public :

    DenseEigenProblem();

    void InitMatrix(Matrix<Tstiff, Prop, Storage>&);
    void InitMatrix(Matrix<Tstiff, Prop, Storage>&, Matrix<Tmass, PropM, StorageM>& );

    void ComputeDiagonalMass();
    void FactorizeCholeskyMass();
    
    void MltCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X);
    void MltCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X);
    
    void SolveCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X);
    void SolveCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X);

    void MltMass(const Vector<T>& X, Vector<T>& Y);    
    
    void MltStiffness(const Vector<T>& X, Vector<T>& Y);
    void MltStiffness(const T& coef_mass, const T& coef_stiff,
		      const Vector<T>& X, Vector<T>& Y);

    void ComputeAndFactoRealMatrix(const Treal&, const Treal& a,
				   const Treal& b, int which);
    
    void ComputeAndFactoRealMatrix(const Tcplx&, const Treal& a,
				   const Treal& b, int which);

    void ComputeAndFactorizeStiffnessMatrix(const Treal& a, const Treal& b,
					    int which_part =
					    EigenProblem_Base<T>::COMPLEX_PART);
    
    void ComputeAndFactorizeStiffnessMatrix(const Tcplx& a, const Tcplx& b,
                                            int which_part =
					    EigenProblem_Base<T>::COMPLEX_PART);
    
    void ComputeSolution(const Vector<Treal>& X, Vector<Treal>& Y);
    void ComputeSolution(const Vector<Tcplx>& X, Vector<Tcplx>& Y);

    void ComputeSolution(const SeldonTranspose& transA,
                         const Vector<Treal>& X, Vector<Treal>& Y);

    void ComputeSolution(const SeldonTranspose& transA,
                         const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    void Clear();
    
  };
  
	
  //! computation of a few eigenvalues for sparse matrices
  template<class T, class MatStiff,
           class MatMass = Matrix<typename ClassComplexType<T>::Treal,
				  Symmetric, ArrayRowSymSparse> >
  class SparseEigenProblem
    : public VirtualEigenProblem<T, typename MatStiff::entry_type,
				 typename MatMass::entry_type>
  {
  protected :
    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;
    
    //! LU factorisation of sparse matrix
    SparseDirectSolver<Treal> mat_lu_real;
    
    //! factorisation of complex system
    SparseDirectSolver<Tcplx> mat_lu_cplx;
    
    //! Cholesky factorisation of mass matrix if required
    SparseCholeskySolver<Treal> chol_facto_mass_matrix;
    
    //! temporary vectors for Cholesky
    Vector<Treal> Xchol_real, Xchol_imag;

    MatMass* Mh;
    MatStiff* Kh;
    
  public :
    
    SparseEigenProblem();
    
    void SelectCholeskySolver(int type);
    
    void InitMatrix(MatStiff&);
    void InitMatrix(MatStiff&, MatMass&);

    void ComputeDiagonalMass();
    void FactorizeCholeskyMass();
    
    template<class Storage>
    void FactorizeCholeskyMass(Matrix<Treal, Symmetric, Storage>& M);
    
    template<class T0, class Prop, class Storage>
    void FactorizeCholeskyMass(Matrix<T0, Prop, Storage>& M);
    
    void MltCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X);
    void MltCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X);
    
    void SolveCholeskyMass(const SeldonTranspose& transA, Vector<Treal>& X);
    void SolveCholeskyMass(const SeldonTranspose& transA, Vector<Tcplx>& X);

    void ComputeAndFactoRealMatrix(const Treal&, const Treal& a,
				   const Treal& b, int which);

    void ComputeAndFactoRealMatrix(const Tcplx&, const Treal& a,
				   const Treal& b, int which);
    
    void ComputeAndFactorizeStiffnessMatrix(const Treal& a, const Treal& b,
					    int which =
					    EigenProblem_Base<T>::COMPLEX_PART);
    
    void ComputeAndFactorizeStiffnessMatrix(const Tcplx& a, const Tcplx& b,
                                            int which =
					    EigenProblem_Base<T>::COMPLEX_PART);
    
    void ComputeSolution(const Vector<Treal>& X, Vector<Treal>& Y);
    void ComputeSolution(const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    void ComputeSolution(const SeldonTranspose& transA,
                         const Vector<Treal>& X, Vector<Treal>& Y);

    void ComputeSolution(const SeldonTranspose& transA,
                         const Vector<Tcplx>& X, Vector<Tcplx>& Y);
    
    void Clear();
    
  };
 
    
  //! list of availables eigenvalue solvers
  class TypeEigenvalueSolver
  {
  public :
    enum {DEFAULT, ARPACK, ANASAZI, FEAST};
    
    static int default_solver;
    
    static int GetDefaultSolver();
    
  };
  
  template<class T, class Prop, class Storage>
  void GetEigenvaluesEigenvectors(EigenProblem_Base<T>& var_eig,
				  Vector<T>& lambda, Vector<T>& lambda_imag,
				  Matrix<T, Prop, Storage>& eigen_vec);
  
}

#define SELDON_FILE_VIRTUAL_EIGENVALUE_SOLVER_HXX
#endif


// Copyright (C) 2009-2010 INRIA
// Author(s): Vivien Mallet
//
// This file is part of the data assimilation library Verdandi.
//
// Verdandi is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Verdandi. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_MODEL_QUADRATICMODEL_HXX


namespace Verdandi
{


    /////////////////////
    // QUADRATIC MODEL //
    /////////////////////


    //! This class is a quadratic model.
    /*! The model is defined as \f$\frac{\mathrm{d}x_i}{\mathrm{d}t} = x^T Q_i
      x + L_i x + b_i\f$, where \f$Q_i\f$ is a matrix, \f$L_i\f$ is the
      \f$i\f$-th row of the matrix \f$L\f$ and \f$b\f$ a vector.
      \tparam T the type of floating-point numbers.
    */
    template <class T>
    class QuadraticModel: public VerdandiBase
    {
    public:
        typedef T value_type;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;
        typedef Matrix<T> tangent_linear_operator;
        typedef Matrix<T> state_error_variance;
        typedef Matrix<T> state_error_variance_reduced;
        typedef Vector<T> state_error_variance_row;
        typedef Vector<T> state;
        typedef Matrix<T> matrix_state_observation;
        typedef Matrix<T> error_variance;
        typedef Vector<T> uncertain_parameter;
        typedef Matrix<T, Symmetric, Seldon::RowSymPacked> parameter_variance;

    protected:

        //! Dimension of the state.
        int Nstate_;

        //! State vector.
        Vector<T> state_;

        //! Source vector.
        state source_;

        //! Should the quadratic term be applied?
        bool with_quadratic_term_;
        //! Should the linear term be applied?
        bool with_linear_term_;
        //! Should the constant term be applied?
        bool with_constant_term_;

        //! Quadratic terms.
        vector<Matrix<T> > S_;

        //! Matrix that defines the linear part of the model.
        Matrix<T> L_;

        //! Vector that defines the constant part of the model.
        Vector<T> b_;

        //! Time step.
        double Delta_t_;

        //! Final time of the simulation.
        double final_time_;

        //! Current time.
        double time_;

        //! Temporary variable that stores S times the state vector.
        Vector<T> S_state_;

        /*** Uncertainty on parameters ***/

        //! Should the quadratic term be perturbed?
        bool is_quadratic_perturbed_;
        //! Should the linear term be perturbed?
        bool is_linear_perturbed_;
        //! Should the constant term be perturbed?
        bool is_constant_perturbed_;

        //! Parameters to be perturbed.
        vector<uncertain_parameter*> parameter_;

        //! List of parameters to be perturbed.
        vector<string> uncertain_parameter_vector_;

        //! Number of parameters to be perturbed.
        int Nparameter_;

        //! Name of the parameters to be perturbed.
        vector<string> parameter_name_;

        //! Correlations between the constant term and the other terms.
        vector<Vector<T> > correlation_;

        //! Name of the probability distribution for the constant term.
        vector<string> pdf_;

        //! Mean of the probability distribution for the constant term.
        vector<Vector<T> > mean_;

        //! Covariance matrix for the constant term.
        vector<Matrix<T, Symmetric, RowSymPacked> > variance_;

        //! PDF parameters for the constant term.
        vector<Vector<T> > optional_parameters_;


        /*** Errors ***/

        //! Variance of the model error.
        error_variance Q_;

        //! Variance of the model error in square root form.
        error_variance Q_sqrt_;

        //! Variance of the state error.
        state_error_variance P_;

        /*! \brief Projector matrix L in the decomposition of the
        background error covariance matrix (\f$B\f$) as a product LUL^T */
        state_error_variance state_error_variance_projector_;

        /*! \brief Reduced matrix U in the decomposition of the
          background error covariance matrix (\f$B\f$) as a product LUL^T */
        state_error_variance_reduced state_error_variance_reduced_;


        //! Variance of the state error in square root form.
        state_error_variance P_sqrt_;

        //! Tangent linear operator (H).
        tangent_linear_operator tangent_linear_operator_;

        //! Index of the row of P currently stored.
        int current_row_;
        //! Value of the row of P currently stored.
        state_error_variance_row state_error_variance_row_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

#ifdef VERDANDI_WITH_MPI
        //! Communicator used inside this model.
        MPI_Comm mpi_communicator_;
#endif

    public:
        // Constructors and destructor.
        QuadraticModel();
        QuadraticModel(string configuration_file);
        ~QuadraticModel();

        // Initializations.
        void Initialize(string configuration_file);
        void InitializeStep();

        // Processing.
        void Forward();
        double ApplyOperator(state& x, bool preserve_state = true);
        double ApplyTangentLinearOperator(state& x);
        tangent_linear_operator& GetTangentLinearOperator();

        bool HasFinished() const;
        void Save();

        void FinalizeStep();
        void Finalize();

        // Access methods.
        T GetDelta_t() const;
        double GetTime() const;
        void SetTime(double time);
        int GetNstate() const;
        int GetNfull_state() const;
        state& GetState();
        void StateUpdated();
        state& GetFullState();
        void FullStateUpdated();
        void SetSource(state& source);

        int GetNparameter();
        uncertain_parameter& GetParameter(int i);
        void ParameterUpdated(int i);
        string GetParameterName(int i);
        Vector<T>& GetParameterCorrelation(int i);
        string GetParameterPDF(int i);
        parameter_variance& GetParameterVariance(int i);
        Vector<T>& GetParameterPDFData(int i);
        string GetParameterOption(int i);

#ifdef VERDANDI_WITH_MPI
        // MPI specific methods.
        void SetMPICommunicator(MPI_Comm);
#endif

        // Errors.
        error_variance& GetErrorVariance();
#ifndef SWIG
        const error_variance& GetErrorVariance() const;
#endif
        error_variance& GetErrorVarianceSqrt();
#ifndef SWIG
        const error_variance& GetErrorVarianceSqrt() const;
#endif
        state_error_variance& GetStateErrorVariance();
#ifndef SWIG
        const state_error_variance& GetStateErrorVariance() const;
#endif
        state_error_variance_row& GetStateErrorVarianceRow(int row);
        state_error_variance& GetStateErrorVarianceSqrt();
#ifndef SWIG
        const state_error_variance& GetStateErrorVarianceSqrt() const;
#endif
        state_error_variance& GetStateErrorVarianceProjector();
        state_error_variance_reduced& GetStateErrorVarianceReduced();

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_QUADRATICMODEL_HXX
#endif

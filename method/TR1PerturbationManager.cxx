// Copyright (C) 2010-2012 INRIA
// Author(s): Kévin Charpentier, Vivien Mallet, Anne Tilloy
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


#ifndef VERDANDI_FILE_METHOD_TR1PERTURBATIONMANAGER_CXX
#define VERDANDI_FILE_METHOD_TR1PERTURBATIONMANAGER_CXX


#include "TR1PerturbationManager.hxx"
#include "BasePerturbationManager.cxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    /*! The seed is initialized from the system clock.
     */
    TR1PerturbationManager
    ::TR1PerturbationManager():
        BasePerturbationManager<TR1PerturbationManager>(), urng_(NULL)
    {
        urng_ = new engine(time(NULL));
        distribution_uniform_ = new distribution_uniform(0.0, 1.0);
        variate_generator_uniform_ =
            new generator_uniform(*urng_, *distribution_uniform_);

        distribution_normal_ = new distribution_normal(0.0, 1.0);
        variate_generator_normal_ =
            new generator_normal(*urng_, *distribution_normal_);
    }


    //! Main constructor.
    /*! Builds the manager and reads option keys in the configuration file.
      \param[in] configuration_file configuration file.
    */
    TR1PerturbationManager
    ::TR1PerturbationManager(string configuration_file):
        BasePerturbationManager<TR1PerturbationManager>(), urng_(NULL)
    {
        urng_ = new engine(time(NULL));
        distribution_uniform_ = new distribution_uniform(0.0, 1.0);
        variate_generator_uniform_ =
            new generator_uniform(*urng_, *distribution_uniform_);
        Initialize(configuration_file);
        distribution_normal_ = new distribution_normal(0.0, 1.0);
        variate_generator_normal_ =
            new generator_normal(*urng_, *distribution_normal_);
    }


    //! Constructor by copy.
    /*! The seed is initialized from the system clock.
     */
    TR1PerturbationManager
    ::TR1PerturbationManager(const TR1PerturbationManager &unUsed):
        BasePerturbationManager<TR1PerturbationManager>(), urng_(NULL)
    {
        urng_ = new engine(time(NULL));
        distribution_uniform_ = new distribution_uniform(0.0, 1.0);
        variate_generator_uniform_ =
            new generator_uniform(*urng_, *distribution_uniform_);

        distribution_normal_ = new distribution_normal(0.0, 1.0);
        variate_generator_normal_ =
            new generator_normal(*urng_, *distribution_normal_);
    }


    //! Destructor.
    TR1PerturbationManager::~TR1PerturbationManager()
    {
        delete variate_generator_uniform_;
        delete variate_generator_normal_;
        delete urng_;
        delete distribution_normal_;
        delete distribution_uniform_;
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the manager.
    /*!
      \param[in] configuration_file configuration file.
    */
    void TR1PerturbationManager::Initialize(string configuration_file)
    {
        VerdandiOps configuration_stream(configuration_file);
        Initialize(configuration_stream);
    }


    //! Initializes the manager.
    /*!
      \param[in] configuration_stream configuration stream.
    */
    void TR1PerturbationManager::Initialize(VerdandiOps&
                                            configuration_stream)
    {
        if(urng_ == NULL)
            urng_ = new tr1::mt19937(time(NULL));
    }


    //! Does nothing.
    void TR1PerturbationManager::Finalize()
    {
    }


    //! Generates a random number with a normal distribution.
    /*!
      \param[in] mean mean of the distribution.
      \param[in] variance variance of the random variable.
      \param[in] parameter vector of parameters. The vector may either be
      empty or contain two clipping parameters \f$ (a, b) \f$. With the
      clipping parameters, for a normal distribution, any random value lies in
      \f$ [\mu - a \sigma, \mu + b \sigma] \f$ where \f$ \mu \f$ is the mean
      of the random variable and \f$ \sigma \f$ is its standard deviation.
      \return A random number following the previously described normal
      distribution.
    */
    double TR1PerturbationManager
    ::Normal(double mean, double variance,
             Vector<double, VectFull>& parameter)
    {

        double value = mean + (*variate_generator_normal_)() * sqrt(variance);
        if (parameter.GetLength() == 2)
            while (value < parameter(0) || value > parameter(1))
                value = mean + (*variate_generator_normal_)()
                    * sqrt(variance);
        else if (parameter.GetLength() != 0)
            throw ErrorArgument("TR1PerturbationManager"
                                "::Normal(double, double, Vector)",
                                "The vector of parameters should be either "
                                "empty or of length 2, but it contains "
                                + to_str(parameter.GetLength())
                                + " element(s).");

        return value;
    }


    //! Generates a random number with a log-normal distribution.
    /*!
      \param[in] mean mean of the distribution.
      \param[in] variance variance of the random variable.
      \param[in] parameter vector of parameters. The vector may either be
      empty or contain two clipping parameters \f$ (a, b) \f$. With the
      clipping parameters, for a normal distribution, any random value lies in
      \f$ [\mu - a \sigma, \mu + b \sigma] \f$ where \f$ \mu \f$ is the mean
      of the random variable and \f$ \sigma \f$ is its standard deviation.
      \return A random number following the previously described normal
      distribution.
    */
    double TR1PerturbationManager
    ::LogNormal(double mean, double variance,
                Vector<double, VectFull>& parameter)
    {
        if (parameter.GetLength() != 0 && parameter.GetLength() != 2)
            throw ErrorArgument("TR1PerturbationManager"
                                "::LogNormal(double, double, Vector)",
                                "The vector of parameters should be either "
                                "empty or of length 2, but it contains "
                                + to_str(parameter.GetLength())
                                + " element(s).");
        return exp(Normal(mean, variance, parameter));
    }


    //! Generates a random number from a uniform distribution.
    /*!
      \param[in] min lower bound of the distribution.
      \param[in] max upper bound of the distribution.
      \return A random number sampled from the uniform distribution U(\a min,
      \a max).
    */
    double TR1PerturbationManager::Uniform(double min, double max)
    {
        double span = (max > min)? max - min: min - max;
        double value = (*variate_generator_uniform_)() * span + min;

        return value;
    }


    //! Generates a random number from a uniform distribution on integers.
    /*!
      \param[in] min lower bound of the distribution.
      \param[in] max upper bound of the distribution.
      \return A random integer sampled from the uniform distribution U(\a min,
      \a max).
    */
    int TR1PerturbationManager::UniformInt(int min, int max)
    {
         std::tr1::uniform_int<int> unif(min, max - 1);
        return unif(*urng_);
    }


    //! Generates a vector of random numbers with a normal distribution.
    /*! Each component of the random vector is generated independently.
      \param[in] variance variance of the random variable.
      \param[in] parameter vector of parameters. The vector may either be
      empty or contain two clipping parameters \f$ (a, b) \f$. With the
      clipping parameters, for a normal distribution, any random value lies in
      \f$ [\mu - a \sigma, \mu + b \sigma] \f$ where \f$ \mu \f$ is the mean
      of the random variable and \f$ \sigma \f$ is its standard deviation.
      \param[out] output the generated random vector.
    */
    template <class T0, class T1,
              class Prop0, class Allocator0>
    void TR1PerturbationManager
    ::Normal(Matrix<T0, Prop0, RowSymPacked, Allocator0> variance,
             Vector<double, VectFull>& parameter,
             Vector<T1, VectFull, Allocator0>& output)
    {
        if (parameter.GetLength() != 0 && parameter.GetLength() != 2)
            throw ErrorArgument("TR1PerturbationManager"
                                "::Normal(string, double, double,"
                                "Vector, Vector)",
                                "The vector of parameters should be either "
                                "empty or of length 2, but it contains "
                                + to_str(parameter.GetLength())
                                + " element(s).");

        Vector<T1, VectFull, Allocator0> sample(output.GetSize());

        int m = variance.GetM();
        Vector<T0, VectFull> diagonal(m);
        for (int i = 0; i < m; i++)
            diagonal(i) = sqrt(variance(i, i));

        GetCholesky(variance);
        Matrix<T1, General, RowMajor> standard_deviation(m, m);
        standard_deviation.Zero();
        for (int i = 0; i < m; i++)
            for (int j = 0; j <= i; j++)
                standard_deviation(i, j) = variance(i, j);

        bool satisfy_constraint = false;

        Vector<T1, VectFull, Allocator0> perturbation(output.GetSize());

        while(!satisfy_constraint)
        {
            perturbation.Zero();
            double value;
            int size = sample.GetSize();
            for (int i = 0; i < size; i++)
            {
                value = (*variate_generator_normal_)();
                if (parameter.GetLength() == 2)
                    while (value < parameter(0) || value > parameter(1))
                        value = (*variate_generator_normal_)();
                sample(i) = value;
            }

            MltAdd(T0(1), standard_deviation, sample,
                   T1(0), perturbation);

            satisfy_constraint = NormalClipping(diagonal, parameter,
                                                perturbation);
        }
        Add(T1(1), perturbation, output);
    }


    //! Generate a random vector with a log-normal distribution.
    /*
      \param[in] variance variance of the log-normal distribution.
      \param[in] parameter vector of parameters. The vector may either be
      empty or contain two clipping parameters \f$ (a, b) \f$. With the
      clipping parameters, for a normal distribution, any random value lies in
      \f$ [\mu - a \sigma, \mu + b \sigma] \f$ where \f$ \mu \f$ is the mean
      of the random variable and \f$ \sigma \f$ is its standard deviation.
      \param[in,out] output output on entry, the mean vector; on exit,
      the sample.
    */
    template <class T0, class Prop0, class Allocator0,
              class T1, class Allocator1>
    void TR1PerturbationManager
    ::LogNormal(Matrix<T0, Prop0, RowSymPacked, Allocator0> variance,
                Vector<double, VectFull>& parameter,
                Vector<T1, VectFull, Allocator1>& output)
    {
        int m = variance.GetM();
        for (int i = 0; i < m; i++)
            output(i) = log(output(i));
        Normal(variance, parameter, output);
        for (int i = 0; i < m; i++)
            output(i) = exp(output(i));
    }


    //! Generate a random vector with a homogeneous normal distribution.
    /*
      \param[in] variance variance of the normal distribution.
      \param[in] parameter vector of parameters. The vector may either be
      empty or contain two clipping parameters \f$ (a, b) \f$. With the
      clipping parameters, for a normal distribution, any random value lies in
      \f$ [\mu - a \sigma, \mu + b \sigma] \f$ where \f$ \mu \f$ is the mean
      of the random variable and \f$ \sigma \f$ is its standard deviation.
      \param[in,out] output output on entry, the mean vector; on exit,
      the sample.
    */
    template <class T0,
              class T1, class Allocator1>
    void TR1PerturbationManager
    ::NormalHomogeneous(T0 variance,
                        Vector<double, VectFull>& parameter,
                        Vector<T1, VectFull, Allocator1>& output)
    {
        T1 value;
        value = Normal(T0(0), variance, parameter);
        for (int i = 0; i < output.GetLength(); i++)
            output(i) += value;
    }


    //! Generates a random vector with a homogeneous log normal distribution.
    /*
      \param[in] variance  variance of the log-normal distribution.
      \param[in] parameter vector of parameters. The vector may either be
      empty or contain two clipping parameters \f$ (a, b) \f$. With the
      clipping parameters, for a log-normal distribution, any random value
      lies in \f$ [\mu - a \sigma, \mu + b \sigma] \f$ where \f$ (\mu, \sigma)
      \f$ are the mean and the standard deviation of the logarithm of the
      random variable.
      \param[out] output output on entry, the median vector; on exit,
      the sample.
    */
    template <class T0,
              class T1, class Allocator1>
    void TR1PerturbationManager
    ::LogNormalHomogeneous(T0 variance,
                           Vector<double, VectFull>& parameter,
                           Vector<T1, VectFull, Allocator1>& output)
    {
        T1 value;
        value = LogNormal(T0(0), variance, parameter);
        for (int i = 0; i < output.GetLength(); i++)
            output(i) *= value;
    }


    //! Tests if a vector satisfies clipping constraints.
    /*!
      \param[in] diagonal diagonal coefficients of the covariance matrix.
      \param[in] permutation vector of permutations done during the Cholesky
      decomposition.
      \param[in] parameter vector of parameters. The vector may either be
      empty or contain two clipping parameters \f$ (a, b) \f$. With the
      clipping parameters, for a normal distribution, any random value lies in
      \f$ [\mu - a \sigma, \mu + b \sigma] \f$ where \f$ \mu \f$ is the mean
      of the random variable and \f$ \sigma \f$ is its standard deviation.
      \param[in] output vector to be tested. This vector was generated using
      a covariance matrix with diagonal \a diagonal.
      \return true if the vector satisfies the constraints.
    */
    template <class T0,
              class T1, class Allocator1>
    bool TR1PerturbationManager
    ::NormalClipping(Vector<T0, VectFull>& diagonal,
                     Vector<double, VectFull>& parameter,
                     Vector<T1, VectFull, Allocator1>& output)
    {
        if (parameter.GetLength() == 0)
            return true;
        if (parameter.GetLength() != 2)
            throw ErrorArgument("TR1PerturbationManager::NormalClipping",
                                "The vector of parameters should be either "
                                "empty or of length 2, but it contains "
                                + to_str(parameter.GetLength())
                                + " element(s).");

        if (diagonal.GetLength() != output.GetLength())
            throw ErrorArgument("TR1PerturbationManager::NormalClipping",
                                "The size of the covariance matrix ("
                                + to_str(diagonal.GetLength())
                                + " x " + to_str(diagonal.GetLength()) + ") "
                                + "is incompatible with that of the output ("
                                + to_str(output.GetLength()) + ").");

        int i = 0;
        T1 value;
        while (i < output.GetLength())
        {
            value = output(i) / diagonal(i);
            if (value < parameter(0) || value > parameter(1))
                return false;
            i++;
        }
        return true;
    }


} // namespace Verdandi.


#endif

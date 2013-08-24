/*
 *  Copyright (c) 2009, 2010, Universitaet Bremen
 *  Copyright (c) 2013 DFKI
 *  All rights reserved.
 *
 *  Author: Christoph Hertzberg <chtz@informatik.uni-bremen.de>
 *  Modified by: Javier Hidalgo Carrio <javier.hidalgo_carrio@dfki.de>
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the Universitaet Bremen nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _MTKWRAP_HPP_
#define _MTKWRAP_HPP_

#include <cassert>/** Assert */
#include <Eigen/Core> /** Eigen */


namespace  localization
{

    /**
     * MtkWrap<M> wraps an MTK-Manifold M to a usckf-compatible manifold.
     * M has to have an enum DOF and implement the methods boxplus and boxminus.
     */
    template<class M>
    struct MtkWrap : public M
    {
            typedef MtkWrap<M> self;
    public:
            //typedef double scalar; // MTK only works with double
            typedef typename M::scalar scalar_type;
            typedef typename M::VectorizedType VectorizedType;

            enum {
                    DOF = M::DOF
            };

            typedef Eigen::Matrix<scalar_type, DOF, 1> vectorized_type;

            MtkWrap(const M &m=M()) : M(m) {}

            /*
             * manifold operator (+)
             *
             */
            self& operator+=(const vectorized_type &delta_state)
            {
                    assert(delta_state.size() == DOF);
                    M::boxplus(delta_state.data());
                    return *this;
            }

            const self operator+(const vectorized_type &delta_state) const
            {
                    self result = *this;
                    result += delta_state;

                    return result;
            }

            /*
             * manifold operator (-)
             */
            const vectorized_type operator-(const self &other) const
            {
                    vectorized_type result;
                    assert(result.rows()==DOF);
                    M::boxminus(result.data(), other);

                    return result;
            }

            bool operator==(const self &other) const
            {
                    vectorized_type diff = (*this) - other;
                    return diff.isZero(1e-12);
            }

            bool operator!=(const self &other) const
            {
                    return !(*this == other);
            }

            vectorized_type getVectorizedState(const VectorizedType type = M::EULER_ANGLES)
            {
                return M::getVectorizedState(type);
            }

    public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

} // namespace localization
#endif /* _MTKWRAP_HPP_ */

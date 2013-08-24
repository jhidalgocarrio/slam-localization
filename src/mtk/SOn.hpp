#ifndef _LOCALIZATION_SON_HPP_
#define _LOCALIZATION_SON_HPP_

#include <mtk/types/vect.hpp>
#include <mtk/src/mtkmath.hpp>

#include <Eigen/Geometry>

namespace localization
{
    namespace MTK
    {
        /**
        * Three-dimensional orientations represented as Quaternion.
        * It is assumed that the internal Quaternion always stays normalized,
        * should this not be the case, call inherited member function @c normalize().
        */
        template<class _scalar = double>
        struct SO3 : public Eigen::Quaternion<_scalar> {
            enum {DOF = 3, DIM = 3};
            typedef _scalar scalar;
            typedef Eigen::Quaternion<scalar> base;
            typedef ::MTK::vect<DIM, scalar> vect_type;


            //! Calculate @c this->inverse() * @c r
            SO3 operator%(const base &r) const {
                    return base::conjugate() * r;
            }

            //! Calculate @c this->inverse() * @c r
            template<class Derived>
            vect_type operator%(const Eigen::MatrixBase<Derived> &vec) const {
                    return base::conjugate() * vec;
            }

            //! Calculate @c this * @c r.conjugate()
            SO3 operator/(const base& r) const {
                    return *this * r.conjugate();
            }

            /**
             * Construct from real part and three imaginary parts.
             * Quaternion is normalized after construction.
             */
            SO3(const scalar& w, const scalar& x, const scalar& y, const scalar& z) : base(w, x, y, z) {
                    base::normalize();
            }

            /**
             * Construct from Eigen::Quaternion.
             * @note Non-normalized input may result result in spurious behavior.
             */
            SO3(const base& src = base::Identity()) : base(src) {}

            /**
             * Construct from rotation matrix.
             * @note Invalid rotation matrices may lead to spurious behavior.
             */
            template<class Derived>
            SO3(const Eigen::MatrixBase<Derived>& matrix) : base(matrix) {}

            //! @name Manifold requirements
            void boxplus(::MTK::vectview<const scalar, DOF> vec, scalar scale=1) {

                Eigen::AngleAxisd angleAxis(Eigen::AngleAxisd(vec[2], Eigen::Vector3d::UnitZ())*
                Eigen::AngleAxisd(vec[1], Eigen::Vector3d::UnitY()) *
                Eigen::AngleAxisd(vec[0], Eigen::Vector3d::UnitX()));

                *this = *this * angleAxis;
            }

            void boxminus(::MTK::vectview<scalar, DOF> res, const SO3<scalar>& other) const {
                    res = SO3::log(*this * other.conjugate());
            }

            friend std::ostream& operator<<(std::ostream &os, const SO3<scalar>& q){
                    return os << q.coeffs().transpose() << " ";
            }
            friend std::istream& operator>>(std::istream &is, SO3<scalar>& q){
                    for(int i=0; i<4; ++i)
                            is >> q.coeffs()[i];
                    q.normalize();
                    return is;
            }

             /**
             * Calculate the inverse of @c exp.
             * Only guarantees that <code>exp(log(x)) == x </code>
             */
            static typename base::Vector3 log(const SO3 &orient){
                    typename base::Vector3 res;
                    ::MTK::log<scalar, 3>(res, orient.w(), orient.vec(), scalar(2), true);
                    return res;
	    }
        };
    }// namespace MTK

}// namespace localization
#endif /* _LOCALIZATION_SON_HPP_ */

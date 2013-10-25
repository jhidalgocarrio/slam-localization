#ifndef _STATE_HPP_
#define _STATE_HPP_

/** MTK library **/
#include <mtk/types/pose.hpp>
#include <mtk/src/SubManifold.hpp>

#include <localization/mtk/SOn.hpp>

#ifndef PARSED_BY_DOXYGEN
//////// internals //////

#define MTK_DEFINE_FEATURE_TYPE(type)\
    typedef std::vector< type > MTK##type##Feature; \
    inline std::ostream& operator<<(std::ostream& __os, const MTK##type##Feature& __var) \
    { \
        __os << "["; \
        for (MTK##type##Feature::const_iterator ii = __var.begin(); ii != __var.end(); ++ii) \
        { \
            __os << " " << *ii; \
        } \
        __os << " ]"; \
        return __os; \
    } \


#define MTK_FEATURE_TYPE(type)\
    MTK##type##Feature \


#endif /* NOT PARSED_BY_DOXYGEN */

namespace localization
{
    // We can't use types having a comma inside AutoConstruct macros :(
    typedef ::MTK::vect<3, double> vec3;
    typedef ::localization::MTK::SO3<double> SO3;
//    typedef std::vector<vec3> featureType;
    MTK_DEFINE_FEATURE_TYPE(vec3)

    struct State
    {
        typedef State self;

        ::MTK::SubManifold<vec3, 0> pos;
        ::MTK::SubManifold<vec3, vec3::DOF + 0> vel;
        ::MTK::SubManifold<SO3, vec3::DOF + vec3::DOF + 0> orient;
        ::MTK::SubManifold<vec3, SO3::DOF + vec3::DOF + vec3::DOF + 0> gbias;
        ::MTK::SubManifold<vec3, vec3::DOF + SO3::DOF + vec3::DOF + vec3::DOF + 0> abias;

        enum
        {
            DOF = vec3::DOF + vec3::DOF + SO3::DOF + vec3::DOF + vec3::DOF + 0
        };

        enum VectorizedType
        {
            EULER_ANGLES = 0,
            ERROR_QUATERNION = 1
        };

        typedef vec3::scalar scalar;
        typedef Eigen::Matrix<scalar, DOF, 1> vectorized_type;

        State ( const vec3& pos = vec3(), const vec3& vel = vec3(), const SO3& orient = SO3(), const vec3& gbias = vec3(), const vec3& abias = vec3() )
            : pos(pos), vel(vel), orient(orient), gbias(gbias), abias(abias)
        {}

        /** @brief set the State from a vectorized type State
         */
        void set (const vectorized_type &vstate, const VectorizedType type = EULER_ANGLES)
        {
            pos = vstate.block<vec3::DOF, 1>(0,0); //! Position
            vel = vstate.block<vec3::DOF, 1>(vec3::DOF,0); //! Velocity
            Eigen::Matrix<scalar, SO3::DOF, 1> orientation =  vstate.block<SO3::DOF, 1>(2*vec3::DOF, 0); //! Orientation

            if (type == EULER_ANGLES)
            {
                orient = Eigen::Quaternion<scalar> (Eigen::AngleAxisd(orientation[2], Eigen::Vector3d::UnitZ())*
                                Eigen::AngleAxisd(orientation[1], Eigen::Vector3d::UnitY()) *
                                Eigen::AngleAxisd(orientation[0], Eigen::Vector3d::UnitX()));
            }
            else
            {
                orient = Eigen::Quaternion<scalar> (1.0, orientation[0], orientation[1], orientation[2]);
            }

            gbias = vstate.block<localization::vec3::DOF, 1>(2*vec3::DOF + SO3::DOF, 0); //! Gyros bias
            abias = vstate.block<localization::vec3::DOF, 1>(3*vec3::DOF + SO3::DOF, 0); //! Acc bias
        }

        int getDOF() const
        {
            return DOF;
        }

        void boxplus(const ::MTK::vectview<const scalar, DOF> & __vec, scalar __scale = 1 )
        {
            pos.boxplus(::MTK::subvector(__vec, &self::pos), __scale);
            vel.boxplus(::MTK::subvector(__vec, &self::vel), __scale);
            orient.boxplus(::MTK::subvector(__vec, &self::orient), __scale);
            gbias.boxplus(::MTK::subvector(__vec, &self::gbias), __scale);
            abias.boxplus(::MTK::subvector(__vec, &self::abias), __scale);
        }

        void boxminus(::MTK::vectview<scalar,DOF> __res, const State& __oth) const
        {
            pos.boxminus(::MTK::subvector(__res, &self::pos), __oth.pos);
            vel.boxminus(::MTK::subvector(__res, &self::vel), __oth.vel);
            orient.boxminus(::MTK::subvector(__res, &self::orient), __oth.orient);
            gbias.boxminus(::MTK::subvector(__res, &self::gbias), __oth.gbias);
            abias.boxminus(::MTK::subvector(__res, &self::abias), __oth.abias);
        }

        friend std::ostream& operator<<(std::ostream& __os, const State& __var)
        {
            return __os << __var.pos << " " << __var.vel << " " << __var.orient << " " << __var.gbias << " " << __var.abias << " " ;
        }

        friend std::istream& operator>>(std::istream& __is, State& __var)
        {
            return __is >> __var.pos >> __var.vel >> __var.orient >> __var.gbias >> __var.abias ;
        }

        /**@brief Create a vectorize state of a single state vector
         * but not in the form of Manifold in the form of error quaternion
         * for the orientation.
         */

        vectorized_type getVectorizedState (const VectorizedType type = EULER_ANGLES)
        {

            State::vectorized_type vstate;

            vstate.block<vec3::DOF, 1>(0,0) = pos; //! Position
            vstate.block<vec3::DOF, 1>(vec3::DOF,0) = vel; //! Velocity
            Eigen::Matrix<scalar, SO3::DOF, 1> orientation; //! Orientation

            if (type == EULER_ANGLES)
            {
                orientation[2] = orient.toRotationMatrix().eulerAngles(2,1,0)[0];//Yaw
                orientation[1] = orient.toRotationMatrix().eulerAngles(2,1,0)[1];//Pitch
                orientation[0] = orient.toRotationMatrix().eulerAngles(2,1,0)[2];//Roll
            }
            else
            {
                orientation << orient.x(), orient.y(), orient.z();
            }

            vstate.block<SO3::DOF, 1>(2*vec3::DOF, 0) = orientation;
            vstate.block<localization::vec3::DOF, 1>(2*vec3::DOF + SO3::DOF, 0) = gbias; //! Gyros bias
            vstate.block<localization::vec3::DOF, 1>(3*vec3::DOF + SO3::DOF, 0) = abias; //! Acc bias

            return vstate;
        }

    };


    struct AugmentedState
    {
        typedef AugmentedState self;

        ::MTK::SubManifold<State, 0> statek; /** Oldest pose state(when first exteroceptive measurement was taken) */
        ::MTK::SubManifold<State, State::DOF + 0> statek_l; /** Pose state (when second exteroceptive measurement was taken) */
        ::MTK::SubManifold<State, State::DOF + State::DOF + 0> statek_i; /** Current Pose state (update to the proprioceptive measurements) */
        ::MTK::SubManifold<MTK_FEATURE_TYPE(vec3), State::DOF + State::DOF + State::DOF + 0> featuresk; /** Features of the measurement took at t=k */
        ::MTK::SubManifold<MTK_FEATURE_TYPE(vec3), State::DOF + State::DOF + State::DOF + 0> featuresk_l; /** Features of the measurement took at t=k+l */

        enum
        {
            DOF = State::DOF + State::DOF + State::DOF + 0
        };

        enum VectorizedType
        {
            EULER_ANGLES = 0,
            ERROR_QUATERNION = 1
        };


        typedef vec3::scalar scalar;
        typedef Eigen::Matrix<scalar, DOF, 1> vectorized_type;

        AugmentedState ( const State& statek = State(),
                const State& statek_l = State(),
                const State& statek_i = State(),
                const MTK_FEATURE_TYPE(vec3)& featuresk = MTK_FEATURE_TYPE(vec3)(),
                const MTK_FEATURE_TYPE(vec3)& featuresk_l = MTK_FEATURE_TYPE(vec3)()
                )
            : statek(statek), statek_l(statek_l), statek_i(statek_i), featuresk(featuresk), featuresk_l(featuresk_l)
            {}

        int getDOF() const
        {
            return DOF;
        }

        int getCurrentDOF() const
        {
            return DOF+featuresk.size()+featuresk_l.size();
        }

        void boxplus(const ::MTK::vectview<const scalar, DOF> & __vec, scalar __scale = 1 )
        {
            statek.boxplus(::MTK::subvector(__vec, &self::statek), __scale);
            statek_l.boxplus(::MTK::subvector(__vec, &self::statek_l), __scale);
            statek_i.boxplus(::MTK::subvector(__vec, &self::statek_i), __scale);

        }

        void boxminus(::MTK::vectview<scalar,DOF> __res, const AugmentedState& __oth) const
        {
            statek.boxminus(::MTK::subvector(__res, &self::statek), __oth.statek);
            statek_l.boxminus(::MTK::subvector(__res, &self::statek_l), __oth.statek_l);
            statek_i.boxminus(::MTK::subvector(__res, &self::statek_i), __oth.statek_i);
        }

        friend std::ostream& operator<<(std::ostream& __os, const AugmentedState& __var)
        {
            __os << "\n" << __var.statek << "\n" << __var.statek_l << "\n" << __var.statek_i << "\n";
            __os << "[";
            for (MTK_FEATURE_TYPE(vec3)::const_iterator ii = __var.featuresk.begin(); ii != __var.featuresk.end(); ++ii)
            {
                __os << " " << *ii;
            }
            __os << " ]\n";
            __os << "[";
            for (MTK_FEATURE_TYPE(vec3)::const_iterator ii = __var.featuresk_l.begin(); ii != __var.featuresk_l.end(); ++ii)
            {
                __os << " " << *ii;
            }
            __os << " ]\n";

            return __os;
        }
        friend std::istream& operator>>(std::istream& __is, AugmentedState& __var)
        {
            return __is >> __var.statek >> __var.statek_l >> __var.statek_i;
        }

        vectorized_type getVectorizedState (const VectorizedType type = EULER_ANGLES)
        {
            AugmentedState::vectorized_type vstate;

            /** statek **/
            vstate.block<State::DOF, 1> (0,0) = statek.getVectorizedState(static_cast<State::VectorizedType>(type));

            /** statek_l **/
            vstate.block<State::DOF, 1> (State::DOF,0) = statek_l.getVectorizedState(static_cast<State::VectorizedType>(type));

            /** statek_i **/
            vstate.block<State::DOF, 1> (2*State::DOF,0) = statek_i.getVectorizedState(static_cast<State::VectorizedType>(type));

            return vstate;
        }

    };

}

#endif /** end of _STATE_HPP_ */

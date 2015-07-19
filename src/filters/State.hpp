#ifndef _STATE_HPP_
#define _STATE_HPP_

/** MTK library **/
#include <mtk/types/pose.hpp>
#include <mtk/src/SubManifold.hpp>

//#include <localization/mtk/SOn.hpp>
#include <mtk/types/SOn.hpp>

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
    typedef ::MTK::SO3<double> SO3;
//    typedef std::vector<vec3> featureType;
//    MTK_DEFINE_FEATURE_TYPE(vec3)

    struct State
    {
        typedef State self;

        ::MTK::SubManifold<vec3, 0> pos;
        ::MTK::SubManifold<SO3, vec3::DOF + 0> orient;
        ::MTK::SubManifold<vec3, vec3::DOF + SO3::DOF + 0> velo;
        ::MTK::SubManifold<vec3, vec3::DOF + SO3::DOF + vec3::DOF + 0> angvelo;

        enum
        {
            DOF = vec3::DOF + SO3::DOF + vec3::DOF + vec3::DOF + 0
        };

        enum VectorizedMode
        {
            EULER_ANGLES = 0,
            ANGLE_AXIS = 1
        };

        typedef vec3::scalar scalar;
        typedef Eigen::Matrix<scalar, DOF, 1> vectorized_type;

        State ( const vec3& pos = vec3(), const SO3& orient = SO3(), const vec3& velo = vec3(), const vec3& angvelo = vec3() )
            : pos(pos), orient(orient), velo(velo), angvelo(angvelo)
        {}

        /** @brief set the State from a vectorized type State
         */
        void set (const vectorized_type &vstate, const VectorizedMode type = ANGLE_AXIS)
        {
            pos = vstate.block<vec3::DOF, 1>(0,0); //! Position
            Eigen::Matrix<scalar, SO3::DOF, 1> axis_angle =  vstate.block<SO3::DOF, 1>(vec3::DOF, 0); //! Orientation

            if (type == EULER_ANGLES)
            {
                orient = Eigen::Quaternion<scalar> (Eigen::AngleAxisd(axis_angle[2], Eigen::Vector3d::UnitZ())*
                                Eigen::AngleAxisd(axis_angle[1], Eigen::Vector3d::UnitY()) *
                                Eigen::AngleAxisd(axis_angle[0], Eigen::Vector3d::UnitX()));
            }
            else
            {
                orient = SO3::exp(axis_angle, 1);
            }

            velo = vstate.block<localization::vec3::DOF, 1>(vec3::DOF + SO3::DOF, 0); //! Linear Velocity
            angvelo = vstate.block<localization::vec3::DOF, 1>(2*vec3::DOF + SO3::DOF, 0); //! Angular Velocity
        }

        void boxplus(const ::MTK::vectview<const scalar, DOF> & __vec, scalar __scale = 1 )
        {
            pos.boxplus(::MTK::subvector(__vec, &self::pos), __scale);
            orient.boxplus(::MTK::subvector(__vec, &self::orient), __scale);
            velo.boxplus(::MTK::subvector(__vec, &self::velo), __scale);
            angvelo.boxplus(::MTK::subvector(__vec, &self::angvelo), __scale);
        }

        void boxminus(::MTK::vectview<scalar,DOF> __res, const State& __oth) const
        {
            pos.boxminus(::MTK::subvector(__res, &self::pos), __oth.pos);
            orient.boxminus(::MTK::subvector(__res, &self::orient), __oth.orient);
            velo.boxminus(::MTK::subvector(__res, &self::velo), __oth.velo);
            angvelo.boxminus(::MTK::subvector(__res, &self::angvelo), __oth.angvelo);
        }

        friend std::ostream& operator<<(std::ostream& __os, const State& __var)
        {
            return __os << __var.pos << " " << " " << __var.orient << " " << __var.velo << " " << __var.angvelo << " " ;
        }

        friend std::istream& operator>>(std::istream& __is, State& __var)
        {
            return __is >> __var.pos >> __var.orient >> __var.velo >> __var.angvelo ;
        }

        /**@brief Create a vectorize state of a single state vector
         */

        vectorized_type getVectorizedState (const VectorizedMode type = ANGLE_AXIS)
        {

            State::vectorized_type vstate;

            vstate.block<vec3::DOF, 1>(0,0) = pos; //! Position
            Eigen::Matrix<scalar, SO3::DOF, 1> orientation; //! Orientation

            if (type == EULER_ANGLES)
            {
                orientation[2] = orient.toRotationMatrix().eulerAngles(2,1,0)[0];//Yaw
                orientation[1] = orient.toRotationMatrix().eulerAngles(2,1,0)[1];//Pitch
                orientation[0] = orient.toRotationMatrix().eulerAngles(2,1,0)[2];//Roll
            }
            else
            {
                orientation << SO3::log(orient);
            }

            vstate.block<SO3::DOF, 1>(vec3::DOF, 0) = orientation;
            vstate.block<localization::vec3::DOF, 1>(vec3::DOF + SO3::DOF, 0) = velo; //! Linear Velocity
            vstate.block<localization::vec3::DOF, 1>(2*vec3::DOF + SO3::DOF, 0) = angvelo; //! Angular Velocity

            return vstate;
        }
    };

    struct SensorState
    {
        typedef SensorState self;

        ::MTK::SubManifold<vec3, 0> pos;
        ::MTK::SubManifold<SO3, vec3::DOF + 0> orient;

        enum
        {
            DOF = vec3::DOF + SO3::DOF + 0
        };

        enum VectorizedMode
        {
            EULER_ANGLES = 0,
            ANGLE_AXIS = 1
        };

        typedef vec3::scalar scalar;
        typedef Eigen::Matrix<scalar, DOF, 1> vectorized_type;

        SensorState ( const vec3& pos = vec3(), const SO3& orient = SO3())
            : pos(pos), orient(orient)
        {}

        /** @brief set the Sensor state from a vectorized type Sensor State
         */
        void set (const vectorized_type &vstate, const VectorizedMode type = ANGLE_AXIS)
        {
            pos = vstate.block<vec3::DOF, 1>(0,0); //! Position
            Eigen::Matrix<scalar, SO3::DOF, 1> axis_angle =  vstate.block<SO3::DOF, 1>(vec3::DOF, 0); //! Orientation

            if (type == EULER_ANGLES)
            {
                orient = Eigen::Quaternion<scalar> (Eigen::AngleAxisd(axis_angle[2], Eigen::Vector3d::UnitZ())*
                                Eigen::AngleAxisd(axis_angle[1], Eigen::Vector3d::UnitY()) *
                                Eigen::AngleAxisd(axis_angle[0], Eigen::Vector3d::UnitX()));
            }
            else
            {
                orient = SO3::exp(axis_angle, 1);
            }
        }

        void boxplus(const ::MTK::vectview<const scalar, DOF> & __vec, scalar __scale = 1 )
        {
            pos.boxplus(::MTK::subvector(__vec, &self::pos), __scale);
            orient.boxplus(::MTK::subvector(__vec, &self::orient), __scale);
        }

        void boxminus(::MTK::vectview<scalar,DOF> __res, const SensorState& __oth) const
        {
            pos.boxminus(::MTK::subvector(__res, &self::pos), __oth.pos);
            orient.boxminus(::MTK::subvector(__res, &self::orient), __oth.orient);
        }

        friend std::ostream& operator<<(std::ostream& __os, const SensorState& __var)
        {
            return __os << __var.pos << " " << " " << __var.orient << " " ;
        }

        friend std::istream& operator>>(std::istream& __is, SensorState& __var)
        {
            return __is >> __var.pos >> __var.orient;
        }

        /**@brief Create a vectorize state of a single state vector
         */

        vectorized_type getVectorizedState (const VectorizedMode type = ANGLE_AXIS)
        {

            SensorState::vectorized_type vstate;

            vstate.block<vec3::DOF, 1>(0,0) = pos; //! Position
            Eigen::Matrix<scalar, SO3::DOF, 1> orientation; //! Orientation

            if (type == EULER_ANGLES)
            {
                orientation[2] = orient.toRotationMatrix().eulerAngles(2,1,0)[0];//Yaw
                orientation[1] = orient.toRotationMatrix().eulerAngles(2,1,0)[1];//Pitch
                orientation[0] = orient.toRotationMatrix().eulerAngles(2,1,0)[2];//Roll
            }
            else
            {
                orientation << SO3::log(orient);
            }

            vstate.block<SO3::DOF, 1>(vec3::DOF, 0) = orientation;

            return vstate;
        }
    };

    template <class _SensorState>
    struct MultiState
    {
        typedef MultiState self;

        ::MTK::SubManifold<State, 0> statek; /** Current Pose state (update to the proprioceptive measurements) */
        std::vector<_SensorState> sensorsk; /** Multi state sensor poses **/

        enum
        {
            SENSOR_DOF = _SensorState::DOF
        };

        enum
        {
            DOF = State::DOF + 0
        };

        enum VectorizedMode
        {
            EULER_ANGLES = 0,
            ANGLE_AXIS = 1
        };


        typedef vec3::scalar scalar;
        typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> vectorized_type;

        typedef State SingleState;

        MultiState (
                const State& statek = State(),
                const  std::vector<_SensorState> &sensorsk = std::vector<_SensorState>()
                )
            : statek(statek), sensorsk(sensorsk)
            {}

        unsigned int getDOF() const
        {
            return State::DOF + (SENSOR_DOF * sensorsk.size());
        }

        /** @brief set the State from a vectorized type State
         */
        void set (const vectorized_type &vstate, const VectorizedMode type = ANGLE_AXIS)
        {
            /** Set state **/
            Eigen::Matrix<scalar, State::DOF, 1> tmp_vstate;
            tmp_vstate = vstate.block(0, 0, State::DOF, 1);
            statek.set(tmp_vstate, localization::State::VectorizedMode(type));

            /** Set sensor poses states **/
            register size_t sensor_idx = 0;
            for (typename std::vector<_SensorState>::iterator it = sensorsk.begin();
                    it != sensorsk.end(); ++it)
            {
                Eigen::Matrix<scalar, _SensorState::DOF, 1> tmp_vsensor;
                tmp_vsensor = vstate.block(State::DOF +(sensor_idx*_SensorState::DOF), 0, _SensorState::DOF, 1);
                it->set(tmp_vsensor, _SensorState::VectorizedMode(type));
                sensor_idx++;
            }

            return;
        }

        void boxplus(MultiState & __state, scalar __scale = 1 )
        {
            State::vectorized_type vectstate;

            vectstate = __state.statek.getVectorizedState();
            statek.boxplus(vectstate.data(), __scale);

            typename std::vector<_SensorState>::iterator it_this_sensorsk = sensorsk.begin();
            typename std::vector<_SensorState>::iterator it_state_sensorsk = __state.sensorsk.begin();
            for ( ;it_state_sensorsk != __state.sensorsk.end(); ++it_state_sensorsk, ++it_this_sensorsk)
            {
                typename _SensorState::vectorized_type vectsensor;
                vectsensor = it_state_sensorsk->getVectorizedState();
                it_this_sensorsk->boxplus(vectsensor.data(), __scale);
            }
        }

        void boxminus(MultiState &__res, const MultiState& __oth) const
        {
            State::vectorized_type vectstate;
            //std::cout<<"in boxminus __res:\n "<<__res<<"\n";
            //std::cout<<"in boxminus __oth:\n "<<__oth<<"\n";
            vectstate = __res.statek.getVectorizedState();
            statek.boxminus(vectstate.data(), __oth.statek);
            __res.statek.set(vectstate);

            typename std::vector<_SensorState>::const_iterator it_this_sensorsk = sensorsk.begin();
            typename std::vector<_SensorState>::iterator it_res_sensorsk = __res.sensorsk.begin();
            typename std::vector<_SensorState>::const_iterator it_oth_sensorsk = __oth.sensorsk.begin();
            for ( ;it_res_sensorsk != __res.sensorsk.end();
                    ++it_res_sensorsk, ++it_this_sensorsk, ++it_oth_sensorsk)
            {
                typename _SensorState::vectorized_type vectsensor;
                vectsensor = it_res_sensorsk->getVectorizedState();
                it_this_sensorsk->boxminus(vectsensor.data(), *it_oth_sensorsk);
                it_res_sensorsk->set(vectsensor);
            }

            //std::cout<<"after in boxminus __res:\n "<<__res<<"\n";
        }

        friend std::ostream& operator<<(std::ostream& __os, const MultiState& __var)
        {
            __os << "\n" << __var.statek << "\n";

            for (typename std::vector<_SensorState>::const_iterator it = __var.sensorsk.begin();
                    it != __var.sensorsk.end(); ++it)
            {
                __os << *it << "\n";
            }

            return __os;
        }

        friend std::istream& operator>>(std::istream& __is, MultiState& __var)
        {
            __is >> __var.statek;

            for (typename std::vector<_SensorState>::iterator it = __var.sensorsk.begin();
                    it != __var.sensorsk.end(); ++it)
            {
                __is >> *it;
            }

            return __is;
        }

        vectorized_type getVectorizedState (const VectorizedMode type = ANGLE_AXIS)
        {
            vectorized_type vstate(static_cast<unsigned int>(DOF), 1);

            /** Statek **/
            vstate.block(0, 0, State::DOF, 1) = statek.getVectorizedState(static_cast<State::VectorizedMode>(type));

            /** Sensors **/
            register size_t sensor_idx = 0;
            for (typename std::vector<_SensorState>::iterator it = sensorsk.begin();
                    it != sensorsk.end(); ++it)
            {
                vstate.block(State::DOF + (sensor_idx*_SensorState::DOF), 0, _SensorState::DOF, 1) = it->getVectorizedState(static_cast<typename _SensorState::VectorizedMode >(type));
                sensor_idx++;
            }

            return vstate;
        }
    };

    template < int _MeasurementDimension >
    struct AugmentedState
    {
        typedef AugmentedState self;

        typedef ::MTK::vect<_MeasurementDimension, double> MeasurementType;

        ::MTK::SubManifold<State, 0> statek; /** Oldest pose state(when first exteroceptive measurement was taken) */
        ::MTK::SubManifold<State, State::DOF + 0> statek_l; /** Pose state (when second exteroceptive measurement was taken) */
        ::MTK::SubManifold<State, State::DOF + State::DOF + 0> statek_i; /** Current Pose state (update to the proprioceptive measurements) */
        MeasurementType featuresk; /** Features of the measurement have taken at t=k */
        MeasurementType featuresk_l; /** Features of the measurement have taken at t=k+l */

        enum
        {
            DOF = State::DOF + State::DOF + State::DOF + 0
        };

        enum VectorizedMode
        {
            EULER_ANGLES = 0,
            ANGLE_AXIS = 1
        };


        typedef vec3::scalar scalar;
        typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> vectorized_type;

        AugmentedState ( const State& statek = State(),
                const State& statek_l = State(),
                const State& statek_i = State(),
                const MeasurementType& featuresk = MeasurementType(),
                const MeasurementType& featuresk_l = MeasurementType()
                )
            : statek(statek), statek_l(statek_l), statek_i(statek_i), featuresk(featuresk), featuresk_l(featuresk_l)
            {}
        /** @brief set the State from a vectorized type State
         */
        void set (const vectorized_type &vstate, const size_t size_featuresk, const size_t size_featuresk_l, const VectorizedMode type = ANGLE_AXIS)
        {
            Eigen::Matrix<scalar, State::DOF, 1> tmp_vstate;
            tmp_vstate = vstate.block(0, 0 ,State::DOF, 1);
            statek.set(tmp_vstate, localization::State::VectorizedMode(type));

            tmp_vstate = vstate.block(State::DOF, 0 ,State::DOF, 1);
            statek_l.set(tmp_vstate, localization::State::VectorizedMode(type));

            tmp_vstate = vstate.block(2*State::DOF, 0 ,State::DOF, 1);
            statek_i.set(tmp_vstate, localization::State::VectorizedMode(type));

            featuresk.resize(size_featuresk, 1);
            featuresk = vstate.block(3*State::DOF, 0, size_featuresk, 1);
            //std::cout<<"set featuresk:\n"<<featuresk<<"\n";

            featuresk_l.resize(size_featuresk_l, 1);
            featuresk_l = vstate.block(3*State::DOF+size_featuresk, 0, size_featuresk_l, 1);
            //std::cout<<"set featuresk_l:\n"<<featuresk_l<<"\n";

            return;
        }

        unsigned int getDOF() const
        {
            return DOF+featuresk.size()+featuresk_l.size();
        }

        void boxplus(AugmentedState & __state, scalar __scale = 1 )
        {
            State::vectorized_type vectstate;

            vectstate = __state.statek.getVectorizedState();
            statek.boxplus(vectstate.data(), __scale);

            vectstate = __state.statek_l.getVectorizedState();
            statek_l.boxplus(vectstate.data(), __scale);

            vectstate = __state.statek_i.getVectorizedState();
            statek_i.boxplus(vectstate.data(), __scale);

            featuresk = featuresk + __state.featuresk;

            featuresk_l = featuresk_l + __state.featuresk_l;
        }

        void boxminus(AugmentedState &__res, const AugmentedState& __oth) const
        {
            State::vectorized_type vectstate;
            //std::cout<<"in boxminus __res:\n "<<__res<<"\n";
            //std::cout<<"in boxminus __oth:\n "<<__oth<<"\n";
            vectstate = __res.statek.getVectorizedState();
            statek.boxminus(vectstate.data(), __oth.statek);
            __res.statek.set(vectstate);

            vectstate = __res.statek_l.getVectorizedState();
            statek_l.boxminus(vectstate.data(), __oth.statek_l);
            __res.statek_l.set(vectstate);

            vectstate = __res.statek_i.getVectorizedState();
            statek_i.boxminus(vectstate.data(), __oth.statek_i);
            __res.statek_i.set(vectstate);

            __res.featuresk = __res.featuresk - __oth.featuresk;
            __res.featuresk_l = __res.featuresk_l - __oth.featuresk_l;

            //std::cout<<"after in boxminus __res:\n "<<__res<<"\n";
        }

        friend std::ostream& operator<<(std::ostream& __os, const AugmentedState& __var)
        {
            __os << "\n" << __var.statek << "\n" << __var.statek_l << "\n"
            << __var.statek_i << "\n" << __var.featuresk.matrix()<<"\n"<< __var.featuresk_l.matrix()<<"\n";

            return __os;
        }
        friend std::istream& operator>>(std::istream& __is, AugmentedState& __var)
        {
            return __is >> __var.statek >> __var.statek_l >> __var.statek_i >> __var.featuresk >> __var.featuresk_l;
        }

        vectorized_type getVectorizedState (const VectorizedMode type = ANGLE_AXIS)
        {
            vectorized_type vstate(this->getDOF(), 1);

            /** statek **/
            vstate.block<State::DOF, 1> (0,0) = statek.getVectorizedState(static_cast<State::VectorizedMode>(type));

            /** statek_l **/
            vstate.block<State::DOF, 1> (State::DOF,0) = statek_l.getVectorizedState(static_cast<State::VectorizedMode>(type));

            /** statek_i **/
            vstate.block<State::DOF, 1> (2*State::DOF,0) = statek_i.getVectorizedState(static_cast<State::VectorizedMode>(type));

            /** featuresk **/
            vstate.block(3*State::DOF, 0, this->featuresk.size(), 1) = featuresk;

            /** featuresk_l **/
            vstate.block(3*State::DOF+this->featuresk.size(), 0, this->featuresk_l.size(), 1) = featuresk_l;

            return vstate;
        }
    };
}

#endif /** end of _STATE_HPP_ */

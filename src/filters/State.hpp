#ifndef _STATE_HPP_
#define _STATE_HPP_

/** MTK library **/
#include <mtk/types/pose.hpp>
#include <mtk/types/SOn.hpp>
#include <mtk/src/SubManifold.hpp>


#ifndef PARSED_BY_DOXYGEN
//////// internals //////

#define MTK_DEFINE_FEATURE_TYPE(type)\
    typedef std::vector< type > MTK##type##Feature; \
    std::ostream& operator<<(std::ostream& __os, const MTK##type##Feature& __var) \
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
    typedef MTK::vect<3, double> vec3;
    typedef MTK::SO3<double> SO3;
    typedef std::vector<vec3> featureType;
    MTK_DEFINE_FEATURE_TYPE(vec3)

    MTK_BUILD_MANIFOLD ( SingleState ,
    (( vec3, pos ))
    (( vec3, vel ))
    (( SO3, orient ))
    (( vec3, gbias ))
    (( vec3, abias ))
    );

    struct VectorState
    {
        typedef VectorState self;

        MTK::SubManifold<SingleState, 0> statek; /** Oldest pose state(when first extereoceptive measurement was taken) */
        MTK::SubManifold<SingleState, SingleState::DOF + 0> statek_l; /** Pose state (when second extereoceptive measurement was taken) */
        MTK::SubManifold<SingleState, SingleState::DOF + SingleState::DOF + 0> statek_i; /** Current Pose state (update to the proprioceptive measurements) */
        MTK::SubManifold<MTK_FEATURE_TYPE(vec3), SingleState::DOF + SingleState::DOF + SingleState::DOF + 0> featuresk; /** Features of the measurement took at t=k */
        MTK::SubManifold<MTK_FEATURE_TYPE(vec3), SingleState::DOF + SingleState::DOF + SingleState::DOF + 0> featuresk_l; /** Features of the measurement took at t=k+l */

        enum
        {
            DOF = SingleState::DOF + SingleState::DOF + SingleState::DOF + 0
        };

        typedef vec3::scalar scalar;

        VectorState ( const SingleState& statek = SingleState(),
                const SingleState& statek_l = SingleState(),
                const SingleState& statek_i = SingleState(),
                const MTK_FEATURE_TYPE(vec3)& featuresk = MTK_FEATURE_TYPE(vec3)(),
                const MTK_FEATURE_TYPE(vec3)& featuresk_l = MTK_FEATURE_TYPE(vec3)()
                )
            : statek(statek), statek_l(statek_l), statek_i(statek_i), featuresk(featuresk), featuresk_l(featuresk_l) {}

        int getDOF() const
        {
            return DOF;
        }

        int getCurrentDOF() const
        {
            return DOF+featuresk.size()+featuresk_l.size();
        }

        void boxplus(const MTK::vectview<const scalar, DOF> & __vec, scalar __scale = 1 )
        {
            statek.boxplus(MTK::subvector(__vec, &self::statek), __scale);
            statek_l.boxplus(MTK::subvector(__vec, &self::statek_l), __scale);
            statek_i.boxplus(MTK::subvector(__vec, &self::statek_i), __scale);

        }

        void boxminus(MTK::vectview<scalar,DOF> __res, const VectorState& __oth) const
        {
            statek.boxminus(MTK::subvector(__res, &self::statek), __oth.statek);
            statek_l.boxminus(MTK::subvector(__res, &self::statek_l), __oth.statek_l);
            statek_i.boxminus(MTK::subvector(__res, &self::statek_i), __oth.statek_i);
        }

        friend std::ostream& operator<<(std::ostream& __os, const VectorState& __var)
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
        friend std::istream& operator>>(std::istream& __is, VectorState& __var)
        {
            return __is >> __var.statek >> __var.statek_l >> __var.statek_i;
        }
    };

}

#endif /** end of _STATE_HPP_ */

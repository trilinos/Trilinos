#include <gtest/gtest.h>
#include <stk_emend/independent_set/IndependentSetFinder.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <vector>

namespace
{

class TestInfo
{
public:
    using ExclusionIdentifierType = size_t;
    using GlobalId = size_t;

    TestInfo(
                            const int owner,
                            const GlobalId priorityAndId,
                            const std::vector<ExclusionIdentifierType> &conflictingIds,
                            const std::vector<int> &procsThatNeedToKnow)
    : mOwner{owner},
      mId{priorityAndId},
      mPriority{priorityAndId},
      mTiePriority{priorityAndId},
      mConflictingGlobalIds{conflictingIds},
      mProcsThatNeedToKnow{procsThatNeedToKnow}
    {}

    TestInfo(
                            const int owner,
                            const GlobalId id,
                            const size_t priority,
                            const size_t tiePriority,
                            const std::vector<ExclusionIdentifierType> &conflictingIds,
                            const std::vector<int> &procsThatNeedToKnow)
    : mOwner{owner},
      mId{id},
      mPriority{priority},
      mTiePriority{tiePriority},
      mConflictingGlobalIds{conflictingIds},
      mProcsThatNeedToKnow{procsThatNeedToKnow}
    {}

    const GlobalId &get_unique_id() const { return mId; }
    int get_owner() const { return mOwner; }
    const std::vector<int> &get_procs_that_need_to_know_about_this_info() const { return mProcsThatNeedToKnow; }
    const std::vector<ExclusionIdentifierType> &get_conflicting_ids() const { return mConflictingGlobalIds; }

    class Comparator
    {
    public:
        bool is_first_higher_priority_than_second(const TestInfo& a,const TestInfo& b) const
        {
            return a.mPriority > b.mPriority;
        }

        bool does_first_win_priority_tie_with_second(const TestInfo& a,const TestInfo& b) const
        {
            return a.mTiePriority > b.mTiePriority;
        }
    };

    friend std::ostream & operator<<(std::ostream & os, const TestInfo &info);

private:
    const int mOwner;
    const GlobalId mId;
    const size_t mPriority;
    const size_t mTiePriority;
    const std::vector<ExclusionIdentifierType> mConflictingGlobalIds;
    const std::vector<int> mProcsThatNeedToKnow;
};

std::ostream & operator<<(std::ostream & os, const TestInfo &info)
{
    os << "Info: mOwner: " << info.mOwner
       << ", mId: " << info.mId
       << ", mPriority: " << info.mPriority
       << ", mTiePriority: " << info.mTiePriority
       << ", mConflictingGlobalIds: {";
    for( const TestInfo::ExclusionIdentifierType &other : info.mConflictingGlobalIds)
        os << other << ", ";
    os << "} mProcsThatNeedToKnow: {";
    for( const int proc : info.mProcsThatNeedToKnow)
        os << proc << ", ";
    os << "}";
    return os;
}

bool operator==(const TestInfo &a, const TestInfo &b)
{
    return a.get_unique_id() == b.get_unique_id() &&
           a.get_owner() == b.get_owner() &&
           a.get_procs_that_need_to_know_about_this_info() == b.get_procs_that_need_to_know_about_this_info() &&
           a.get_conflicting_ids() == b.get_conflicting_ids();
}

TEST(TestInfoOutputOperator, callIt_noUncoveredLines)
{
    const TestInfo info(0, 17, {1}, {0});
    std::ostringstream ss;
    ss << info;
    EXPECT_FALSE(ss.str().empty());
}

class TestIndependentSetFinder : public ::testing::Test
{
protected:
    TestInfo::Comparator mComparator;
    const stk::ParallelMachine mComm {MPI_COMM_WORLD};
    std::vector<int> mAllProcs;

    TestIndependentSetFinder()
    {
        for (int i{0}; i < stk::parallel_machine_size(mComm); ++i)
            mAllProcs.push_back(i);
    }

    void test_independent_set(const std::vector<TestInfo> &infos, const std::vector<TestInfo> &goldIndependentInfos)
    {
        std::vector<TestInfo> independentInfos = independent_set::IndependentSetFinder<TestInfo>::find_independent_set(
                                infos,
                                mComparator,
                                mComm);
        EXPECT_EQ(goldIndependentInfos, independentInfos);
    }

    TestInfo make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(
                            const TestInfo::GlobalId id,
                            const size_t priority,
                            const size_t tiePriority,
                            const std::vector<TestInfo::ExclusionIdentifierType> &exclusionIds)
    {
        return TestInfo(0, id, priority, tiePriority, exclusionIds, mAllProcs);
    }

    TestInfo make_info_owned_on_p0_that_all_procs_need_to_know(
                            const TestInfo::GlobalId id,
                            const std::vector<TestInfo::ExclusionIdentifierType> &exclusionIds)
    {
        return TestInfo(0, id, exclusionIds, mAllProcs);
    }
};

TEST_F(TestIndependentSetFinder,noInfos_findIndependentSet_emptyIndependentSet)
{
    test_independent_set({}, {});
}

TEST_F(TestIndependentSetFinder,oneInfo_findIndependentSet_1InfoIndependentSet)
{
    std::vector<TestInfo> infos{ make_info_owned_on_p0_that_all_procs_need_to_know(17, {}) };
    test_independent_set(infos, infos);
}

TEST_F(TestIndependentSetFinder,twoConflictingInfos_findIndependentSet_1InfoIndependentSet)
{
    const std::vector<TestInfo> infos
    {
        make_info_owned_on_p0_that_all_procs_need_to_know(17, {13}),
        make_info_owned_on_p0_that_all_procs_need_to_know(18, {13}),
    };
    test_independent_set(infos, {infos[1]});
}

TEST_F(TestIndependentSetFinder,twoNonConflictingInfos_findIndependentSet_2InfoIndependentSet)
{
    const std::vector<TestInfo> infos
    {
        make_info_owned_on_p0_that_all_procs_need_to_know(19, {}),
        make_info_owned_on_p0_that_all_procs_need_to_know(20, {}),
    };
    test_independent_set(infos, infos);
}

TEST_F(TestIndependentSetFinder,threeConflictingInfos_findIndependentSet_1InfoIndependentSet)
{
    const std::vector<TestInfo> infos
    {
        make_info_owned_on_p0_that_all_procs_need_to_know(21, {1}),
        make_info_owned_on_p0_that_all_procs_need_to_know(22, {1,2}),
        make_info_owned_on_p0_that_all_procs_need_to_know(23, {2}),
    };
    test_independent_set(infos, {infos[2]});
}

TEST_F(TestIndependentSetFinder,threeConflictingInfosBut2LowestPrioritiesAreTied_findIndependentSet_2InfoIndependentSetIncludingTheHighestPriorityAndTheLowestThatDoesntConflictWithHighest)
{
    const std::vector<TestInfo> infos
    {
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(21, 1, 8, {1}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(22, 1, 9, {1,2}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(23, 2, 0, {2}),
    };
    test_independent_set(infos, {infos[0], infos[2]});
}

TEST_F(TestIndependentSetFinder,twoTiedInfosThatDontConflictWithHigherPriorityInfosAnd2ndWinsTie_findIndependentSet_choose2ndWithTieBreaker)
{
    const std::vector<TestInfo> infos
    {
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(20, 1, 7, {1}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(21, 1, 8, {1}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(22, 1, 9, {1,2}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(23, 2, 0, {2}),
    };
    test_independent_set(infos, {infos[1], infos[3]});
}

TEST_F(TestIndependentSetFinder,twoTiedInfosThatDontConflictWithHigherPriorityInfosAnd1stWinsTie_findIndependentSet_choose1stWithTieBreaker)
{
    const std::vector<TestInfo> infos
    {
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(20, 1, 8, {1}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(21, 1, 7, {1}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(22, 1, 9, {1,2}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(23, 2, 0, {2}),
    };
    test_independent_set(infos, {infos[0], infos[3]});
}

TEST_F(TestIndependentSetFinder,FourChainedConflictingInfosWithBottom3TiedAndAlternatingTieBreaker_find_independentSet_ChooseTieBreaker)
{
    const std::vector<TestInfo> infos
    {
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(20, 99, 15, {1}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(21, 99, 10, {1,2}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(22, 99, 20, {2,3}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(23, 100, 0, {3}),
    };
    test_independent_set(infos, {infos[0], infos[3]});
}

TEST_F(TestIndependentSetFinder,SplitInfoChainWhereOneBranchIsHigherAndOtherLowerPriorityThanJoint_findIndependentSet_bothBranchesChosen)
{
    const std::vector<TestInfo> infos
    {
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(19, 99, 5, {1}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(20, 99, 15, {2}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(21, 99, 10, {1,2,3}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(22, 99, 20, {3,4}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(23, 100, 0, {4})
    };
    test_independent_set(infos, {infos[0], infos[1], infos[4]});
}

TEST_F(TestIndependentSetFinder,SplitInfoChainWhereOneBranchIsHigherPriorityAndOtherChainOnOtherProcIsLowerPriorityThanJoint_findIndependentSet_bothBranchesChosen)
{
    if(stk::parallel_machine_size(mComm) >= 2)
    {
        const std::vector<TestInfo> infos
        {
            TestInfo(1, 19, 99, 5, {1}, mAllProcs),
            make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(20, 99, 15, {2}),
            make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(21, 99, 10, {1,2,3}),
            make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(22, 99, 20, {3,4}),
            make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(23, 100, 0, {4})
        };
        test_independent_set(infos, {infos[0], infos[1], infos[4]});
    }
}

TEST_F(TestIndependentSetFinder,TwoChainsThatMeetAtLowestPriorityTiedInfo_findIndependentSet_lowestInfoInAfterResolvingBothChains)
{
    const std::vector<TestInfo> infos
    {
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know( 19, 100, 5, {1}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(119, 100, 5, {101}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know( 20,  99, 5, {1,2}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(120,  99, 5, {101,102}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know(  1,  99, 1, {2,102}),
    };
    test_independent_set(infos, {infos[0], infos[1], infos[4]});
}

TEST_F(TestIndependentSetFinder,threePartiallyConflictingInfos_findIndependentSet_2InfoIndependentSet)
{
    const std::vector<TestInfo> infos
    {
        make_info_owned_on_p0_that_all_procs_need_to_know(21, {1,2}),
        make_info_owned_on_p0_that_all_procs_need_to_know(22, {1}),
        make_info_owned_on_p0_that_all_procs_need_to_know(23, {2}),
    };
    test_independent_set(infos, {infos[1], infos[2]});
}

TEST_F(TestIndependentSetFinder,twoConflictingInfosOnDifferentProcs_findIndependentSet_1InfoIndependentSet)
{
    if(stk::parallel_machine_size(mComm) >= 2)
    {
        const std::vector<TestInfo> infos
        {
            {0, 17, {13}, mAllProcs},
            {1, 18, {13}, mAllProcs},
        };
        test_independent_set(infos, {infos[1]});
    }
}

TEST_F(TestIndependentSetFinder,firstInfoIsBetterThanOneNeighborAndLosesTieToAnOutNeighbor_findIndependentSet_firstInfoShouldBeIn)
{
    const std::vector<TestInfo> infos
    {
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know( 1, 10, 1, {1,2}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know( 2, 9, 2, {1}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know( 3, 10, 2, {2,3}),
        make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know( 4, 11, 0, {3}),
    };
    test_independent_set(infos, {infos[0], infos[3]});
}

class TestIsBetterOrWinsTieInIndependentSet : public TestIndependentSetFinder
{
protected:
    void test_is_better_or_wins_tie_with_unknown_neighbor(
                            bool isInfoTiedWithNeighbor,
                            bool doesInfoWinTieOverNeighbor,
                            bool isNeighborUnknown,
                            bool goldIsBetter)
    {
        const size_t infoPriority = 10;
        const size_t neighborInfoPriority = isInfoTiedWithNeighbor ? infoPriority : infoPriority-1;
        const size_t infoTiedPriority = 10;
        const size_t neighborInfoTiedPriority = doesInfoWinTieOverNeighbor ? infoTiedPriority-1 : infoTiedPriority+1;
        const independent_set::IndependentSetFinder<TestInfo>::InOutUnknownStatus neighborStatus =
                                isNeighborUnknown ?
                                independent_set::IndependentSetFinder<TestInfo>::UNKNOWN_STATUS :
                                independent_set::IndependentSetFinder<TestInfo>::OUT;
        const std::vector<TestInfo> infos
        {
            make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know( 1, infoPriority, infoTiedPriority, {1}),
            make_info_with_id_and_priortity_owned_on_p0_that_all_procs_need_to_know( 2, neighborInfoPriority, neighborInfoTiedPriority, {1}),
        };

        const std::vector<independent_set::IndependentSetFinder<TestInfo>::InOutUnknownStatus> infoInOutStatus
        {
            independent_set::IndependentSetFinder<TestInfo>::UNKNOWN_STATUS,
            neighborStatus
        };

        const independent_set::IndependentSetInfoToInfoGraph<TestInfo> conflictingInfoGraph{infos};
        bool isBetter = independent_set::IndependentSetFinder<TestInfo>::is_better_than_neighbors_or_wins_tie_breaker_with_unknown_neighbors(
                                0,
                                infos,
                                mComparator,
                                conflictingInfoGraph,
                                infoInOutStatus);
        EXPECT_EQ(goldIsBetter, isBetter);
    }

    const bool infoIsTiedWithNeighbor{true};
    const bool infoIsBetterThanNeighbor{false};
    const bool infoWinsTieOverNeighbor{true};
    const bool infoLosesTieOverNeighbor{false};
    const bool neighborIsUnknown{true};
    const bool neighborIsOut{false};
};

TEST_F(TestIsBetterOrWinsTieInIndependentSet,infoThatIsBetterThanNeighborButWouldLoseTieBreakerIfChecked_isBetterOrWinsTie_true)
{
    test_is_better_or_wins_tie_with_unknown_neighbor(infoIsBetterThanNeighbor, infoLosesTieOverNeighbor, neighborIsUnknown, true);
}

TEST_F(TestIsBetterOrWinsTieInIndependentSet,infoThatIsTiedAndWouldLoseTieBreakerWithUnknownNeighbor_isBetterOrWinsTie_false)
{
    test_is_better_or_wins_tie_with_unknown_neighbor(infoIsTiedWithNeighbor, infoLosesTieOverNeighbor, neighborIsUnknown, false);
}

TEST_F(TestIsBetterOrWinsTieInIndependentSet,infoThatIsTiedAndWouldWinTieBreakerWithUnknownNeighbor_isBetterOrWinsTie_true)
{
    test_is_better_or_wins_tie_with_unknown_neighbor(infoIsTiedWithNeighbor, infoWinsTieOverNeighbor, neighborIsUnknown, true);
}

TEST_F(TestIsBetterOrWinsTieInIndependentSet,infoThatIsTiedAndWouldLoseTieBreakerWithOutNeighbor_isBetterOrWinsTie_true)
{
    test_is_better_or_wins_tie_with_unknown_neighbor(infoIsTiedWithNeighbor, infoLosesTieOverNeighbor, neighborIsOut, true);
}
}

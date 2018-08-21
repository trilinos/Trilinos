#include <gtest/gtest.h>

#include <stk_util/util/RemoveIntersection.hpp>


TEST(RemoveIntersection, disjoint)
{
    std::vector<int> v1 = {1, 3, 2, 4};
    std::vector<int> v2 = {5, 6, 7, 8};

    stk::util::remove_intersection(v1, v2);

    std::vector<int> v1gold = {1, 2, 3, 4};
    std::vector<int> v2gold = {5, 6, 7, 8};
    EXPECT_EQ(v1gold, v1);
    EXPECT_EQ(v2gold, v2);
}

TEST(RemoveIntersection, someIntersection)
{
    std::vector<int> v1 = {1, 2, 4, 3, 5};
    std::vector<int> v2 = {5, 6, 7, 8};

    stk::util::remove_intersection(v1, v2);

    std::vector<int> v1gold = {1, 2, 3, 4};
    std::vector<int> v2gold = {6, 7, 8};
    EXPECT_EQ(v1gold, v1);
    EXPECT_EQ(v2gold, v2);
}

TEST(RemoveIntersection, sameVectors)
{
    std::vector<int> v1 = {5, 2, 3, 4, 1};
    std::vector<int> v2 = {1, 2, 3, 4, 5};

    stk::util::remove_intersection(v1, v2);

    std::vector<int> v1gold;
    std::vector<int> v2gold;
    EXPECT_EQ(v1gold, v1);
    EXPECT_EQ(v2gold, v2);
}

struct foo {
    int ival;
    double dval;
};

bool operator==(const foo& lhs, const foo& rhs) {
    return lhs.ival == rhs.ival;
}
bool operator==(int lhs, const foo& rhs) {
    return lhs == rhs.ival;
}

struct foocomp {
    bool operator()(int lhs, const foo& rhs) const {
        return lhs < rhs.ival;
    }
    bool operator()(const foo& lhs, int rhs) const {
        return lhs.ival < rhs;
    }
    bool operator()(const foo& lhs, const foo& rhs) const {
        return lhs.ival < rhs.ival;
    }
    bool operator()(int lhs, int rhs) const {
        return lhs < rhs;
    }
};

TEST(RemoveIntersection, someCrazyIntersection)
{
    std::vector<foo> v1 = { {1, 1.0}, {2, 2.0}, {4, 4.0}, {3, 3.0}, {5, 5.0}};
    std::vector<int> v2 = {5, 6, 7, 8};

    stk::util::remove_intersection(v1, v2, foocomp());

    std::vector<foo> v1gold = {{1, 1.0}, {2, 2.0}, {3, 3.0}, {4, 4.0}};
    std::vector<int> v2gold = {6, 7, 8};
    EXPECT_EQ(v1gold, v1);
    EXPECT_EQ(v2gold, v2);
}


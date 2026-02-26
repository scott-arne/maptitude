#include <gtest/gtest.h>
#include "maptitude/Residue.h"

#include <map>
#include <unordered_set>

using namespace Maptitude;

TEST(ResidueTest, Construction) {
    Residue r("ALA", 123, "A", " ");
    EXPECT_EQ(r.name, "ALA");
    EXPECT_EQ(r.number, 123);
    EXPECT_EQ(r.chain, "A");
    EXPECT_EQ(r.insert_code, " ");
}

TEST(ResidueTest, DefaultInsertCode) {
    Residue r("GLY", 45, "B");
    EXPECT_EQ(r.insert_code, " ");
}

TEST(ResidueTest, ToString) {
    Residue r("ALA", 123, "A", " ");
    std::string s = r.ToString();
    EXPECT_NE(s.find("ALA"), std::string::npos);
    EXPECT_NE(s.find("123"), std::string::npos);
    EXPECT_NE(s.find("A"), std::string::npos);
}

TEST(ResidueTest, Equality) {
    Residue a("ALA", 123, "A", " ");
    Residue b("ALA", 123, "A", " ");
    Residue c("GLY", 123, "A", " ");

    EXPECT_EQ(a, b);
    EXPECT_NE(a, c);
}

TEST(ResidueTest, Ordering) {
    Residue a("ALA", 1, "A");
    Residue b("GLY", 2, "A");
    Residue c("ALA", 1, "B");

    // Same chain, different number
    EXPECT_LT(a, b);
    // Different chain
    EXPECT_LT(a, c);
}

TEST(ResidueTest, UseAsMapKey) {
    std::map<Residue, double> scores;
    Residue r("ALA", 123, "A");
    scores[r] = 0.95;
    EXPECT_EQ(scores.size(), 1u);
    EXPECT_NEAR(scores[r], 0.95, 1e-10);
}

TEST(ResidueTest, Hash) {
    std::unordered_set<Residue> residues;
    residues.insert(Residue("ALA", 1, "A"));
    residues.insert(Residue("GLY", 2, "A"));
    residues.insert(Residue("ALA", 1, "A"));  // Duplicate

    EXPECT_EQ(residues.size(), 2u);
}

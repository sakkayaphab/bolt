#include <gtest/gtest.h>
#include "caller/filemanager.h"
#include "caller/evidence.h"
#include "caller/smithwaterman.h"
#include <bits/stdc++.h>
#include "caller/alignment.h"
#include "caller/readdepthanalysis.h"
#include "bwt/bwt.h"
#include "caller/stringsearch.h"
#include "caller/stringsearchalignment.h"

TEST(FilePathTest, getFilePathName)
{
   std::string samplepath = "file.bam";
   std::string refpath = "file.fasta";
   std::string output = "tempout";
   FileManager fp(samplepath, refpath, output);

   ASSERT_EQ("file.bam", fp.getSamplePath());
   ASSERT_EQ("file.fasta", fp.getReferencePath());
   ASSERT_EQ("tempout", fp.getOutputPath());
}

// TEST(TestSplitRead, testIndelInReads)
// {
//    std::string samplepath = "/data/users/duangdao/kan/sample/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.chr1_231869234-231869675.bam";
//    std::string refpath = "file.fasta";
//    std::string output = "tempout";
//    FileManager fp(samplepath, refpath, output);

//    ASSERT_EQ("file.bam", fp.getSamplePath());
//    ASSERT_EQ("file.fasta", fp.getReferencePath());
//    ASSERT_EQ("tempout", fp.getOutputPath());
// }

//
//TEST(FilePathTest, getFilePathNameGen)
//{
//    std::string samplepath = "file.bam";
//    std::string refpath = "file.fasta";
//    std::string output = "tempout";
//    FileManager fp(samplepath, refpath, output);
//
//    ASSERT_EQ(output + "/analysis/evidence", fp.getEvidencePath());
//}
//
//TEST(Evidence, convertEvidence)
//{
//    Evidence e;
//    e.setEvidenceByString("chr1\t6941739\t.\t.\t.\t.\t.\tEND=6951962;SVTYPE=DEL;SVLEN=10378;CIPOS=-167,167;CIEND=-155,155;BOLT_MQL=60,60,60,60;CHR2=chr2;");
//    ASSERT_EQ(6941739, e.getPos());
//    ASSERT_EQ(6951962, e.getEnd());
//    ASSERT_EQ("chr1", e.getChr());
//    ASSERT_EQ("chr2", e.getEndChr());
//    ASSERT_EQ(4, e.getFrequency());
//    ASSERT_EQ("DEL", e.getVariantType());
//}

// TEST(SmithWaterman, testRealignmnet)
// {
//     // std::string ref = "TCTCACTCTGCTACCCGGGATGGAGTGCAGTGGTGCAATCTCAGCTCACTGTAACCCCTGCCTCCTGGATTCAAGCAATTCTCCTGTCTCAGCCTCCCGAGTAGCTGGGGTTACAGGTGCCATCAGACCCGGCTAATTTTTGTATTTTTAGTAGAGACAGGGTTCACCATGTTGGCCAAGCTGGTCTTGAACTCTTGACCTCAAGTGATCCACCCGCCTCAGTCTCCCAAAGTGCTGGGATTACAGGTGTGAGCCACCACACCCAGCAGTATTAGTTGTTTTATATTGTAGACTCTTAAATGATGTGGACTATGTGTTGTAAAATCCATACCCTGAAGTAAAATGCATTTAGATGGCCCCCACTTCCCTTTTTAACTTAAATCCTTTTTCCCCCTTGCAGTCAGTGCATTTTTATATTTGTACATACTTTGTGAGAAATGTGTAAAGGAAATATTTTTGCATCTAATTGTTCTAACTTCCAAAGATTCTTTACTGAGTTAAAAAAACTTTTTTTTTAGTCTAAACTGAATTTTGTCTGTCTTTCCCCTACTCAGGCCCAGATTGCCTTCCTACAGGGAGAAAGGAAGGGCCAAGAAAATTTGAAGAAGGATCTTGTGAGGAGGATCAAAATGTTGGAGTATGCTCTTAAACAGAAAAGGTAATTCAGTAAAATGAAAAGTGATATTCTTTTTTTGTTTGTTTGTTTTGAGACAGAGTTTCATTCTTGTTGCCCAGGCTGGGGTGAAATGGTGTGATCTTGGCTCACTGCAACCTCCACCTCCCTGTTTCAAGCAATTCTCCTGCCTCAGCCTCCCCAGTAGCTGGGATTATAGGCATGCACCACCACACCCAGCTAATTTTGTATTTTTAGTAGAGACGGGGTTTCTCCATGTTGGTCAGGCTGGTCTCGAACTCCTGACTTCAGGTGATCTGCCCACCTCGGCCTCCTAAAGTGCTGGCATTACAGG";
//     // std::string seq = "GGACAAAGGAAGGGCCAAGAAAATTTGAAGAAGGATCTTGTGAGGAGGATCAAAATG";
//     // std::string ref = "GACCCTCTCACTTCTCCACTGAGTGCCCCGTCTCACTATCCTGCACCCACCCCATCTCCTCTATCCAAACAGGGGGCCCAGTAAGCCTGCCCACTACAGCAAGAAGGAACAGGGCGAGTGCATCCTAAAGGATCCCATCCACGTAACAGTACATGGACAAGAGCAGAAAAAGGAAAGGCAAAATAAAGGTCCCTTGTCATCACGAAAACTAAGATTCTAGACAGGAGTGTTCCCAGGATGGCTGGTCCTCCAGCTCCTCCTGAACACTTTCAGCAGCTCTGTGGTTGCTGGGGACATGACGGTCTGCTGGGGCTCCACAGGGTGACTGGCTGGAGACTGGAAGATTTAGTCTCAGTTTTCCTTCCTATACCAACACCCTGGACCTCCTCCTGCTTTTTTTCCCCCAGAGGCTCCACCTTTAACACAAAGAGGCAAATTCCAGAGAAACATGCTGCCAGCACAGCTGGGCTGGGAGAGGCCCCATCCTCAGGTACCTGACCAAAGCCCCAGCCGCCACAGAGGAGGGCTCGGCCCCCTCTCCTGCCCTGCAAAAGGTGCTATCATCTTTGGTTTTCTTCTTTTCCTCTCTAGCTTCCTCTCTTCCTGCCTCTACCCCAAGGTTTACTGAGGTACAATCAATAAACAAAAAAATTTATATATTCAAGGTGTAAAATGTGACGTTTTAATATATAGAGAGTTGACCCTTGAAGAAGGCAAGAGTTGGGGCACTGATACCCTGCATAGTCAAAAATCCTCATACGACTTTTTAGTCCCACAAAACTTAACTATTATAATAGCCTACTGTTGACCAAAGCCTTCCTGATAACATAAACAGTTGATTAACAGGTATTTTGTGTTATATGTATTATATACTCTATTCTTAAGATAAGCTAGAGAAAAGAAAATGTTATTACAATCATAGGGAAGAGAAAATAGATTTACTGTTCATTAAGTGAAAGTGGATCATCATAAAGGTCTTCATCTTCATTGTCTTTACGTTGAATAGGCTGAAGAGGAGGAGGAAGGGAAAGGGTTGGTCTTGCTGTCTCA";

//     // std::string ref = "CTCGCAGC";
//     // std::string seq = "XXXCTCGCOG";
//     // std::string seq = "XXXCTCGCAG";
//     // std::string seq = "COCGCAGXXX";
//     // std::string seq = "OTCGCAGXXX";
//     // std::string ref = "XXXXXXXXXXTTAAGTATATAAAAAAAAGAAAAAAAGAAAAAAAATCTCTCACCTTTTACCAAAGTGTTGGGATTATGGGTGTGAGCAACCACATCTGGTCCTTTTTTTTT";
//     // std::string ref = "TCTCACTCTGCTACCCGGGATGGAGTGCAGTGGTGCAATCTCAGCTCACTGTAACCCCTGCCTCCTGGATTCAAGCAATTCTCCTGTCTCAGCCTCCCGAGTAGCTGGGGTTACAGGTGCCATCAGACCCGGCTAATTTTTGTATTTTTAGTAGAGACAGGGTTCACCATGTTGGCCAAGCTGGTCTTGAACTCTTGACCTCAAGTGATCCACCCGCCTCAGTCTCCCAAAGTGCTGGGATTACAGGTGTGAGCCACCACACCCAGCAGTATTAGTTGTTTTATATTGTAGACTCTTAAATGATGTGGACTATGTGTTGTAAAATCCATACCCTGAAGTAAAATGCATTTAGATGGCCCCCACTTCCCTTTTTAACTTAAATCCTTTTTCCCCCTTGCAGTCAGTGCATTTTTATATTTGTACATACTTTGTGAGAAATGTGTAAAGGAAATATTTTTGCATCTAATTGTTCTAACTTCCAAAGATTCTTTACTGAGTTAAAAAAACTTTTTTTTTAGTCTAAACTGAATTTTGTCTGTCTTTCCCCTACTCAGGCCCAGATTGCCTTCCTACAGGGAGAAAGGAAGGGCCAAGAAAATTTGAAGAAGGATCTTGTGAGGAGGATCAAAATGTTGGAGTATGCTCTTAAACAGAAAAGGTAATTCAGTAAAATGAAAAGTGATATTCTTTTTTTGTTTGTTTGTTTTGAGACAGAGTTTCATTCTTGTTGCCCAGGCTGGGGTGAAATGGTGTGATCTTGGCTCACTGCAACCTCCACCTCCCTGTTTCAAGCAATTCTCCTGCCTCAGCCTCCCCAGTAGCTGGGATTATAGGCATGCACCACCACACCCAGCTAATTTTGTATTTTTAGTAGAGACGGGGTTTCTCCATGTTGGTCAGGCTGGTCTCGAACTCCTGACTTCAGGTGATCTGCCCACCTCGGCCTCCTAAAGTGCTGGCATTACAGG";
//     // std::string ref = "XXXXXXXXXXCAGAGGTTGCAGTGAGTGGAGAAAGGGAGCGCAGACGCGAG";
//     // std::string ref = "GAGAGAGAGAGA";

//     // std::string seq = "XXXXXTTAAGTATATAA";

//     // std::string seq = "OOOOOTTAAGTATATAA";
//     // SmithWaterman reAlign(&ref, 10, false);
//     // std::vector<SmithWaterman::ScoreAlignment> result = reAlign.findEndToStart(seq);

//     // std::string seq = "TTAAGTATATAAOOOOO";
//     // std::string seq = "CTCGCAGC";

//     // std::string seq = "CAGAGGTTGCAGTGAGTGGAGXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
//     //  std::string seq = "CAGAGGTTGCAGTGAGTGGAGTTAGCAACCCTGAAAGACCAGCAAGTAGAGCTTCGTCCCTCCCTTTTCATGAACCTTAATGTTGTAAGGGAGGCCAAGGAGTGGCCTGGTGAGACTGCACAGCC";
//     // std::string seq = "GAGAAAGGAAGGGCCAAGAAAATTTGAAGAAGGATCTTGTGAGGAGGATCAAAATG";

//     // std::cout << reAlign.getRefPos() << std::endl;
//     // std::vector<SmithWaterman::ScoreAlignment> result = reAlign.findStartToEnd(seq,10);
//     // // std::vector<SmithWaterman::ScoreAlignment> result = reAlign.findStartToEnd(seq,10);

//     // start to end
//     std::string ref = "ATACGGGAGGCTGAGGCAGGAGAATTGCTTGAATCTGGGAGGCAGAGGTTGCAGTGAGTGGAGATCGCACCACTGCACTCCAGCCTGGGAGACAGTGTGAGACTCTGTCTCAAAAAAATAAAAGTAAAATAAAATAAAATAATCACTCTGAATGTAGTGCAGAGAATCTTCCAGGGCAGGGCCAGCAGACCACAGCCCATGGCAGGTTCTGGCCCACTGCCTGTTTCTATATAGCTAGCAAGCTAAGAATGCTTTGTATGTTTTTAAGTGGTTGAAAAAAAATCAGTATTTTATTTCACATGAAAATTATATGAAATTCAAATTTCAGTGTCCAGTTCCAGTTCTATTGAAACAGAGCCACACTCTTGTGTGTATTATCTGAGGCCGCTTTTGAGCTGCAACAGCCAAGCTAAGT";
//     std::string seq = "ATACGGGAGGCTGAGGCAGGAGAATTGCTTGAATCTTGGAGGCAGAGGTTGCAGTGAGTGTAGTTAGCAACCCAGAAAGACCAGCAAGGAGAGCTTCGTCCCTCCCTTTTCATGAACCTGAATGT";
//     // std::string seq = "CCATCTTGTXXXXXXX";
//     // std::string ref = "TTCAAACCCATCTTGTOOOOOOOO";
//     SmithWaterman reAlign(&ref, 0, true);
//     std::vector<SmithWaterman::ScoreAlignment> result = reAlign.findStartToEnd(seq);

//     std::cout << "----------" << std::endl;
//     for (auto n : result)
//     {
//         std::cout << std::endl;
//         std::cout << "pos : " << n.pos << " end : " << n.end
//                   << " posseq : " << n.posseq
//                   << " endseq : " << n.endseq
//                   << " lastscore : " << n.lastscore
//                   << " scorepattern : " << n.scorepattern << std::endl;
//     }

//     // end to start
//     // std::string seq = "XXXXXXCCATCTTGT";
//     // std::string ref = "TTCAAACCCATCTTGT";
//     // SmithWaterman reAlign(&ref, 10, false);
//     // std::vector<SmithWaterman::ScoreAlignment> result = reAlign.findEndToStart(seq);

//     // std::cout << "----------" << std::endl;
//     // for (auto n : result)
//     // {

//     //     std::cout << std::endl;
//     //     std::cout << "pos : " << n.pos << " end : " << n.end
//     //               << " posseq : " << n.posseq
//     //               << " endseq : " << n.endseq
//     //               << " lastscore : " << n.lastscore
//     //               << " scorepattern : " << n.scorepattern << std::endl;
//     // }
// }

// TEST(Aligment, alignDeletionTargetAtStart)
// {
//     std::string ref = "AGGGTGCTAGGTCCTGGGGATACAGCAGTGAGGTCCTTGGTCTCCCAGGGCTTGCTCGGGAACCAGGTGTCCCTCCCATGGGATTCTCGTGTGCTCCCTGGTTATAGCAAGTGCTGTGCTGTGTTTTCGTGATCTGGCTACATGTCTGTTCTCCCCACTAGACCAAGGAGCTTCTCAAGGAGAGAGTCTGAGTCTTCCATTTCTGTATCCCATAACACCTAGTGTTGGGTATATGGAAGGTTCTTTAGAACTGAATAAATGAACTAAAGGGGAGAACAGACCCAGGCCTGCTGATCCCAGGATCAATATGAAATGGGGCAAAGGAGTATTGAGGCAGCTTTTCAGATTCAAAAGCCAAGCTAGCAACAAGTCCCTGGTACAGGGTCTGTGGCTACTGTCAAGGACTGGGCTGTGTGGCCTGGGGACCAACTCACTCCTCTTTTCCTGCCAGTGTGTAGGAGCGGATCCAGGGGTTGGGCACAGACAGCCTGGGGGCCAGG";
//     std::string seq = "CTGGACTACACCTAAAACCAAAACCACTGAGTAAGATTTTTTCATTTTGTGAAAGTTCTGCCAATTTTGCTTAAGTAGAACTGAATAAATGAACTAAAGG";

//     StringSearchAlignment ssa;
//     ssa.setReference(ref);
//     ssa.setSVType("DELEND");
//     ssa.setPosReference(66537691);
//     ssa.buildReference();
//     StringSearchConfig ssc;
//     std::vector<StringSearch::Score> result = ssa.alignDeletionTargetAtStart(&seq,&ssc);
//     for (auto n : result)
//     {
//         int32_t mPos = n.posseq + 66530331;
//         int32_t mEnd = n.pos;

//         std::cout << mPos << " = " << mEnd << std::endl;
//     }
// }

// TEST(Aligment, alignDeletionTargetAtEnd)
// {
//     std::string ref = "OOOOOGGAGAGAGGGXXXXXXXX";
//     std::string seq = "GGAGAGAGGGOOOO";

//     std::cout << "ref : " << ref << std::endl;
//     std::cout << "seq : " << seq << std::endl;

//     StringSearchAlignment ssa;
//     ssa.setReference(ref);
//     ssa.setSVType("DELSTART");
//     ssa.setPosReference(10);
//     ssa.buildReference();
//     ssa.alignDeletionTargetAtEnd(&seq);
// }

// TEST(Aligment, testDuplicationStart)
// {
//     std::string ref = "ACAAAGGAAGGGAGAGAGGG";
//     std::string seq = "ACAAAGGAAGGOOOOOOOOOO";

//     std::cout << "ref : " << ref << std::endl;
//     std::cout << "seq : " << seq << std::endl;

//     StringSearchAlignment ssa;
//     ssa.setReference(ref);
//     ssa.setSVType("DUPSTART");
//     ssa.setPosReference(0);
//     ssa.buildReference();
//     ssa.alignDuplicationTargetAtStart(&seq);
// }

// TEST(Aligment, testDuplicationEnd)
// {
//     std::string ref = "AACGGAGAGAGGG";
//     std::string seq = "OOOOOGGAGAGAGGG";

//     std::cout << "ref : " << ref << std::endl;
//     std::cout << "seq : " << seq << std::endl;

//     StringSearchAlignment ssa;
//     ssa.setReference(ref);
//     ssa.setSVType("DUPEND");
//     ssa.setPosReference(0);
//     ssa.buildReference();
//     ssa.alignDuplicationTargetAtEnd(&seq);
// }

// TEST(Aligment, alignTranslocationTargetAtStartSCS)
// {
//     std::string ref = "ACAAAGGAAGGGAGAGAGGG";
//     std::string seq = "ACAAAGGAAGGOOOOOOOOOO";

//     std::cout << "ref : " << ref << std::endl;
//     std::cout << "seq : " << seq << std::endl;

//     StringSearchAlignment ssa;
//     ssa.setReference(ref);
//     ssa.setSVType("TRASTART");
//     ssa.setPosReference(0);
//     ssa.buildReference();
//     ssa.alignTranslocationTargetAtStartSCS(&seq);
// }

// TEST(Aligment, alignTranslocationTargetAtStartSCE)
// {
//     std::string ref = "AGGGAGAGAGGG";
//     std::string seq = "OOOGGAGAGAGGG";

//     std::cout << "ref : " << ref << std::endl;
//     std::cout << "seq : " << seq << std::endl;

//         StringSearchAlignment ssa;
//     ssa.setReference(ref);
//     ssa.setSVType("TRASTART");
//     ssa.setPosReference(0);
//     ssa.buildReference();
//     ssa.alignTranslocationTargetAtStartSCE(&seq);
// }

// TEST(Aligment, alignTranslocationTargetAtEndSCS)
// {
//     std::string ref = "ACAAAGGAAGGGAGAGAGGG";
//     std::string seq = "ACAAAGGAAGGOOOOOOOOOO";

//     std::cout << "ref : " << ref << std::endl;
//     std::cout << "seq : " << seq << std::endl;

//     StringSearchAlignment ssa;
//     ssa.setReference(ref);
//     ssa.setSVType("TRAEND");
//     ssa.setPosReference(0);
//     ssa.buildReference();
//     ssa.alignTranslocationTargetAtEndSCS(&seq);
// }

// TEST(Aligment, alignTranslocationTargetAtEndSCE)
// {
//     std::string ref = "AGGGAGAGAGGG";
//     std::string seq = "OOOGGAGAGAGGG";

//     std::cout << "ref : " << ref << std::endl;
//     std::cout << "seq : " << seq << std::endl;

//         StringSearchAlignment ssa;
//     ssa.setReference(ref);
//     ssa.setSVType("TRAEND");
//     ssa.setPosReference(10);
//     ssa.buildReference();
//     ssa.alignTranslocationTargetAtEndSCE( &seq);
// }

TEST(Aligment, alignInversionTargetAtStartSCS)
{
   std::string ref = "ATTATTATAGAAAGTTCTTTTTATCAAACCAAACTCTACCTCCCAATAATTTCTGTCCACTGATCTCAGTTCTGATCACTGGGACTTTATAAAATTAGATCAAATCCTGTATCTATTTGGCACCCATGTAATATTAGAAGAGGGTTATATTCTGAGTCGTCTCTTCTGCAAACCAAAACTGTAATCTGCTTCAAATCTCCTTCATAAAAAAAGGTCGGCCAGGTGTGGTGGCTCATGCCTATAATGCAATGCCAGCACTTTAGGAGGCTGAGGTGGGCGGACTGCTTGAGCCCAAGAGTTCAAGACCAGCCTGGACAATATGGCAAAACTCCGTCTCTACAAAAAATACAAAAAAATTAGACCAGCACAGTGGTGCATGCCTGTAGTCCCAGCTACTTGGGAGGCTGGGGCGGGAGGATCCTTGAGCTCAGGAGGCAGAGGTTGCAGTGAGCCGAGATTGTGCTGCTGCACTCCAGCCTGGGTGACAGAGTGAGACCTTATCTCCAAAAAATAATAATAAGACAAGGTTTCTGGACAATGACTTCTATCACTCAACTCAGGATGTACCCTAGTTTAAAAATCTCTTCTTAAACATGACACCCAGAATACTCTCGCAGGTGTACAGTAACCATTGCCAGGTGCCTGGGACTACTGCTTCCACTTTCTGAACTTGAGATCAAAAAAGCCCCCACTGAATTGCTGTTGAGCCAGTTTTCCCACACCCTATCCTTATACAACTGGTTATTTTTCATCTAAATGCAGGACTTTACAGGTACTGCTTTTAAGTGTTTTTGTTTTTTGTTTTTTCAGATAGGGTCTCACTCTGTCACTCAGGCTGGAATGCAGTGGCATGATCTTGGCTCACTGCAACCTCTGCCTCCTGGGTTCAAGTGATTCTCGTGCCTCAACCTCCTGAGTAGCTGGGACTACAGGAGCACACCACCATGCCTGCCTAAATTTTGTATTTTTAGTAGAGGCGGGGTTTCACCATGTTGGCCAGGCTGTTCTTGAACTCCTGACTTCAGGTAGATCTGCCCGCCTTGGCCTTCCAAAGTGCTGGGATTATAGGTGTGAGCCATCGCGCCCGGCCTTTGTCTCTATTTTTGAACATTAGAATTACATTATCTCTCATATAATACATGCTAATATACCAATTAGAACAGAAGTTTTGATTGCTAGACTTGATTTTGTTTACAATTTGGTACAGATCAGGTGTGTATGTTTACTTTTTATGGTTTTTCTTTCATTCCCCTTGCTGATAGCAAGTAGCAGCTACTGTCAG";
   std::string seq = "AAAAATTACCTGGGCGTAGTAGTGTGAGCCTGTAATCCCAGCTACTCAGGAGGCTGAGGCATGAGAATCGCTTGAACCCAGGAGGCAGAGGTTGCAGTGAG";

   std::cout << "ref : " << ref << std::endl;
   std::cout << "seq : " << seq << std::endl;

   StringSearchAlignment ssa;
   ssa.setReference(ref);
   ssa.setSVType("INVSTART");
   ssa.setPosReference(12545664);
   ssa.buildReference();

   StringSearchConfig ssc;
   ssc.setAllowMissMatch(2);
   ssc.setMaxContinueMissMatch(1);
   ssc.setMaxAllowAlign(28);

   std::vector<StringSearch::Score> result = ssa.alignInversionTargetAtStartSCS(&seq, &ssc);
}

// TEST(Aligment, alignInversionTargetAtStartSCE)
// {
//     std::string ref = "TTTTTGCAGGGGGATTTTT";
//     std::string seq = "GGGGGTCCCCCTGC";

//     std::cout << "ref : " << ref << std::endl;
//     std::cout << "seq : " << seq << std::endl;

//     StringSearchAlignment ssa;
//     ssa.setReference(ref);
//     ssa.setSVType("INVSTART");
//     ssa.setPosReference(10);
//     ssa.buildReference();
//     ssa.alignInversionTargetAtStartSCE( &seq);
// }

// TEST(Aligment, alignInversionTargetAtEndSCS)
// {
//     std::string ref = "OOOOOATGGCAGGGGGAC";
//     std::string seq = "GTCCCCCTGCCAOOOOO";

//     std::cout << "ref : " << ref << std::endl;
//     std::cout << "seq : " << seq << std::endl;

//     StringSearchAlignment ssa;
//     ssa.setReference(ref);
//     ssa.setSVType("INVEND");
//     ssa.setPosReference(10);
//     ssa.buildReference();
//     ssa.alignInversionTargetAtStartSCS( &seq);
// }

// TEST(Aligment, testDeletionShortReadStart)
// {
//     std::string ref = "ATAGATCCCTTACCTGGTTTCTTGTGCAATCCCCACTTTTATACTGTTGAGTTTCCCTGATGAGTTTTCATATATATAATCTCTGTGAACTTGGATCAAATTTGGCTTAAATGCCTATTTACCACTTTGCTAAGACAGAAATGTTTGCGCTCAGTCCCTGAATTCCATTGTAGATGTGTTGGCATTTGCTTCCAAACCATCCTTCTGGTTTTAGCTTGCCCTGTGTCTGTGGTCATTTGACATTAGAATAGCCTTCTGTGATCTTCCGGGATCTTAACACCTAAATATTTTCCCTGCTTGTACAGCTAATGAGTTGTGCAGCCTAAGGAAACATGGGGTAGTCAGTTTTCATCAGTATTTCAACTAGGCGTTTTTGCCATCTCATAGCTACTTCTCAAGCAATTACAATTTAAGTCAAGGTAGTAGCACCTGCATCTTCAAAAAACAAAACAAAAAGTTTTTTATCTTCAGAAATAAGGAAAGAAGCTTATTGGACTTCCCAATAGAAAGTCCACCAGCCAATTTTAATACATAAATAATCTATGATTAGCATTAGTATATGAAAAAAATAGTAGCAGCTCAGATTTTATGTACTTCTGCAGTTTAAGAAAGTGATGCCACAGACAACATCTCATTTGGTTTCTTAGCAACCCTGAAAGACCAGCAAGGAGAGCTTCGTCCCTCCCTTTTCATGAACCTGAATGTTGTAAGGGAGGCCAAGGAGTGGCCTGGTGAGACTGCACAGCCAGAGTTAGAGCCCAGGTCCACGTGAGTCATTGTACAGACACATCCTGAGCACCAGCCCACCCACTTCTGCCATCAGGAGGCTTTGAAACTGCTGTGGTTCCCCTCTTAAAATGGCTAACCTCTACTTCAGCAACCTTATTTATGCATTCATATACATAATCATTCTTTAAAGTCTAGTAGACAATTTTCACATTTTGTTTCAGCCCAGATCTTTACTAAATTTTGTGAGATCAGCATTGTTGGGCTTGAATTAACCAAACTTCGCCATACTGTGCCTTCAGAAATCCCAGATTGGGCGGGTTG";
//     std::string seq = "GCTATACGGGAGGCTGAGGCAGGAGAATTGCTTGAATCTGGGAGGCAGAGGTTGCAGTGAGTGGAGTTAGCAACCCTGAAAGACCAGCAAGGAGAGCTTCGTCCCTCCCTTTTCAAGAACCTGAA";

//     std::cout << "ref : " << ref << std::endl;
//     std::cout << "seq : " << seq << std::endl;

//     Alignment alignment(&ref);
//     alignment.genarateMatrix("DELEND");
//     // alignment.setPosReference(10);
//     alignment.alignDeletionTargetAtStart(&seq);
// }

// TEST(Aligment, alignDeletionTargetAtStart)
// {
//     std::string ref = "GGAGAAGAAGAGAGG";
//     std::string seq = "OOOOAGAAGAAGAGAGG";

//     std::cout << "ref : " << ref << std::endl;
//     std::cout << "seq : " << seq << std::endl;
//     // n.pos : 2 n.end : 15 n.posseq : 4 n.endseq : 17 n.pattern : MMMMMMMMMMMMM
//     Alignment alignment(&ref);
//     alignment.genarateMatrix("DELSTART");
//     // alignment.setPosReference(10);
//     alignment.alignDeletionTargetAtStart(&seq);
// }

// TEST(Aligment, alignDeletionTargetAtStart)
// {
//     std::string ref = "ATAGATCCCTTACCTGGTTTOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO";
//     std::string seq = "ATAGATCCCTTACCTGGTTT";
//                 // seq = "MMMMMMMMMMMMMMMMMMMMMM"
//     //  n.pos : 2 n.end : 15 n.posseq : 0 n.endseq : 13 n.pattern : MMMMMMMMMMMMM

//     std::cout << "ref : " << ref << std::endl;
//     std::cout << "seq : " << seq << std::endl;

//     Alignment alignment(&ref);
//     alignment.genarateMatrix("DELSTART");
//     alignment.setPosReference(0);
//     alignment.alignDeletionTargetAtStart(&seq);
// }

// TEST(Aligment, alignDeletionTargetAtEnd)
// {
//     std::string ref = "GAATCTGGGAGGCAGAGGTTGCAGTGAGTGGAGATCGCACCACTGCACTCCAGCCTGGGAGACAGTGTGAGACTCTGTCTCAAAAAAATAAAAGTAAAATAAAATAAAATAATCACTCTGAATGTAGTGCAGAGAATCTTCCAGGGCAGGGCCAGCAGACCACAGCCCATGGCAGGTTCTGGCCCACTGCCTGTTTCTATATAGCTAGCAAGCTAAGAATGCTTTGTATGTTTTTAAGTGGTTGAAAAAAAATCAGTATTTTATTTCACATGAAAATTATATGAAATTCAAATTTCAGTGTCCAGTTCCAGTTCTATTGAAACAGAGCCACACTCTTGTGTGTATTATCTGAGGCCGCTTTTGAGCTGCAACAGCCAAGCTAAGT";
//     std::string seq = "AATCTGGGAGGCAGAGGTGGCAGTGAGTGGAGTTAGCAACCCTGAAAGACCAGCAAGGAGAGCTTCGTCCCTCCCTTTTCATGAACCTGAATGTTGTAAGGGAGGCCAAGGAGTGGCCTGGTGAG";

//     std::cout << "ref : " << ref << std::endl;
//     std::cout << "seq : " << seq << std::endl;

//     Alignment alignment(&ref);
//     alignment.genarateMatrix("DELEND");
//     alignment.setPosReference(10);
//     alignment.alignDeletionTargetAtEnd(&seq);
// }

// TEST(Aligment, readdepthAnalysis)
// {
//     const char *ref = "helloworld";
//     StringSearch ss;
//     ss.setReference(ref);
//     ss.buildHashTable();
//     std::string seq = "hello";
//     std::vector<StringSearch::Score> scoreresult = ss.searchStartToEnd(&seq);
//     for (auto m : scoreresult)
//     {
//         std::cout
//         << "pos : " << m.pos << "\t"
//         << "end : " << m.end << "\t"
//         << "posseq : " << m.posseq << "\t"
//         << "endseq : " << m.endseq << "\t"
//         << std::endl;
//     }
//     // std::cout << std::endl;
// }

// TEST(Aligment, readdepthAnalysis2)
// {
//     const char *ref = "helloworld";
//     StringSearch ss;
//     ss.setReference(ref);
//     ss.buildHashTable();
//     std::string seq = "world";
//     std::vector<StringSearch::Score> scoreresult = ss.searchEndToStart(&seq);
//     for (auto m : scoreresult)
//     {
//         std::cout
//         << "pos : " << m.pos << "\t"
//         << "end : " << m.end << "\t"
//         << "posseq : " << m.posseq << "\t"
//         << "endseq : " << m.endseq << "\t"
//         << std::endl;
//     }
//     std::cout << std::endl;
// }
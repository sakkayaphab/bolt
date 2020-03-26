#include "refineinsertion.h"
#include <stdlib.h>
#include "smithwaterman.h"
#include <graph/dynamicgraph.h>
#include <chrono>

RefineInsertion::RefineInsertion() {
    variantresult.setVariantType("INS");
}

void RefineInsertion::execute() {
    std::cout << evidence.convertToVcfString() << std::endl;

    std::chrono::steady_clock sc;   // create an object of `steady_clock` class
    auto start = sc.now();     // start timer

    variantresult.setChr(evidence.getChr());
    variantresult.setEndChr(evidence.getEndChr());
    prepareBamReader();

    if (evidence.getMark() == "SINS") {
        variantresult = evidence;
        variantresult.setQuailtyPass(true);
        return;
    }

    if (evidence.getMark() == "SR") {
        variantresult = evidence;
        variantresult.setQuailtyPass(true);
        return;
    }

//    std::cout << "END" << std::endl;

    first();

    auto end = sc.now();       // end timer (starting & ending is done by measuring the time at the moment the process started & ended respectively)
    auto time_span = static_cast<std::chrono::duration<double>>(end - start);   // measure time span between start & end
    if (time_span.count()>=10) {
        std::cout<<"Time : "<< time_span.count() << std::endl;
    }
}

void RefineInsertion::first() {
    std::string findRange = convertRangeToString(evidence.getChr(), evidence.getPos() + evidence.getCiPosLeft(),
                                                 evidence.getPos() + evidence.getCiPosRight());

    if (evidence.getPos() + evidence.getCiPosRight() - evidence.getPos() - evidence.getCiPosLeft() < 150) {
        return;
    }

    const char *range = findRange.c_str();
    refineStartToEnd(range);

    RefineInsertion::convertMapSC();
    RefineInsertion::clearMapSC();

//    if (variantresult.isQuailtyPass() == false && evidence.getMark() != "MATEUNMAPPED") {
//        evidence.setMark("UNMERGE");
//        RefineInsertion::findBreakpoint();
//        RefineInsertion::filterBreakpoint();
//    } else {
        RefineInsertion::findBreakpoint();
        RefineInsertion::filterBreakpoint();
//    }
}

void RefineInsertion::refineStartToEnd(const char *range) {
    hts_itr_t *iter = NULL;

    iter = sam_itr_querys(bam_index, bam_header, range);
    if (iter == NULL)
        return;
    read = bam_init1();
    readparser.setBamHeader(bam_header);
    readparser.setBamRead(read);

    int32_t startDist = evidence.getPos() + evidence.getCiPosRight();
    int32_t endDist = evidence.getPos() + evidence.getCiPosLeft();

    while (sam_itr_next(inFile, iter, read) >= 0) {

        if (readparser.isUnmapped()) {
            if (readparser.getSequence().find("N") != std::string::npos) {
                continue;
            }

            vUnMappedRead.push_back(readparser.getSequence());
            continue;
        }

        if (readparser.isNotPassingFilters()) {
            continue;
        }

        if (readparser.isPCR()) {
            continue;
        }

        if (readparser.isSupplementaryAlignment()) {
            continue;
        }

        // 10000
        if (readparser.getPos() < startDist - 10000) {
            break;
        }

        // 20000
        if (readparser.getPos() > endDist + 10000) {
            break;
        }

        auto cigar = readparser.getCigar();
        if (cigar.size() <= 1) {
            continue;
        }

        if (readparser.getSequence().find("N") != std::string::npos) {
            continue;
        }

        if (cigar.at(cigar.size() - 1).getOperatorName() == 'S' && cigar.at(cigar.size() - 1).getLength() >= 4) {
            mapSCEnd[readparser.getEnd()].addMapQ(readparser.getMapQuality());
            mapSCEnd[readparser.getEnd()].addLongMapping(cigar.at(cigar.size() - 1).getLength());
            mapSCEnd[readparser.getEnd()].setPosition(readparser.getEnd());
            mapSCEnd[readparser.getEnd()].addSeqList(
                    readparser.getEdgeSeqFromEndSeq(cigar.at(cigar.size() - 1).getLength()));
        }

        if (cigar.at(0).getOperatorName() == 'S' && cigar.at(0).getLength() >= 4) {
            mapSCStart[readparser.getPos()].addMapQ(readparser.getMapQuality());
            mapSCStart[readparser.getPos()].addLongMapping(cigar.at(0).getLength());
            mapSCStart[readparser.getPos()].setPosition(readparser.getPos());
            mapSCStart[readparser.getPos()].addSeqList(readparser.getEdgeSeqFromStartSeq(cigar.at(0).getLength()));
        }
    }

    hts_itr_destroy(iter);

    return;
}

void RefineInsertion::filterBreakpoint() {
    // std::cout << "# filterBreakpoint : " << vectorBP.size() << std::endl;
    std::sort(vectorBP.begin(), vectorBP.end());

    int maxFreq = 0;
    int maxLongMatch = 0;
    int maxscore = 0;

    for (BreakpointPosition n : vectorBP) {

        int32_t averagePos = 0;
        if (n.pos > n.end) {
            averagePos = n.pos;
        } else {
            averagePos = (n.end + n.pos) / 2;
        }

        if (n.frequency <= 2) {
            continue;
        }

        if (evidence.getMark() == "") {
            if (n.longmatch<50) {
                continue;
            }

//            if (evidence.getFrequency() <= 3) {
//                continue;
//            }
//
//            if (getMaxUInt8FromVector(n.mappingqualitylist) < 10) {
//                continue;
//            }
//
//            if (getMaxUInt8FromVector(n.mappingqualitylist) < 60) {
//                if (n.longmatch <= getDivider(samplestat->getReadLength(), 15, 100, 15)) {
//                    continue;
//                }
//                // continue;
//            } else {
//                if (n.longmatch <= getDivider(samplestat->getReadLength(), 20, 100, 15)) {
//                    continue;
//                }
//                // continue;
//            }
        } else if (evidence.getMark() == "MATEUNMAPPED") {
            if (evidence.getFrequency() <= 4) {
                continue;
            }

            if (n.frequency <= 4) {
                continue;
            }

            if (getMaxUInt8FromVector(n.mappingqualitylist) < 20) {
                continue;
            }

            if (evidence.getMaxMapQ() < 40) {
                continue;
            }

            if (n.longmapstart <= getDivider(samplestat->getReadLength(), 30, 100, 15)) {
                continue;
            }

            if (n.longmapend <= getDivider(samplestat->getReadLength(), 30, 100, 15)) {
                continue;
            }
        }
//
//        if (n.longmatch < getDivider(samplestat->getReadLength(), 15, 100, 15)) {
//            continue;
//        }

        int score = (n.frequency) * (2 * n.longmatch);

        if (score <= maxscore) {
            continue;
        }

        maxscore = score;

        variantresult.setPos(averagePos);
        variantresult.setEnd(averagePos);
        variantresult.setFrequency(n.frequency);
        variantresult.setRPMapQ(*evidence.getMapQVector());
        variantresult.LNGMATCH = n.longmatch;

        variantresult.setMapQList(n.mappingqualitylist);
        variantresult.setChr(evidence.getChr());
        variantresult.setEndChr(evidence.getEndChr());
        variantresult.setQuailtyPass(true);
        variantresult.setMark(evidence.getMark());
    }
}

void RefineInsertion::findBreakpoint() {

//    std::chrono::steady_clock sc;   // create an object of `steady_clock` class
//    auto start = sc.now();     // start timer


    int kmer = 10;
    std::vector<CountRefineSeq> mergeStart;
    std::vector<CountRefineSeq> mergeEnd;

    std::cout << " SCS : " << vectorSCStart.size()
              << " SCE : " << vectorSCEnd.size()
              << " UNMAP : " << vUnMappedRead.size()
              << std::endl;


    for (InsertionPositionDetail n : vectorSCStart) {
        mergeStart.clear();
        if (n.getLongMapping() < kmer) {
            continue;
        }

        bool added;


        mergeStart = mergeString(n, true);



        for (InsertionPositionDetail m : vectorSCEnd) {
            mergeEnd.clear();

            if (m.getLongMapping() < kmer) {
                continue;
            }

            mergeEnd = mergeString(m, true);

            if (m.getPosition()-n.getPosition()>50) {
                continue;
            }

            if (checkBetween(n.getPosition(), m.getPosition(),-samplestat->getReadLength(),samplestat->getReadLength())) {
                if (m.getPosition() - n.getPosition() > 200) {
                    std::cout << n.getPosition() << " " << m.getPosition() << " = " << m.getPosition() - n.getPosition()
                              << std::endl;
                }
            } else {
                continue;
            }


            DynamicGraph dynGraph;
            dynGraph.setKmer(kmer);

            for (int s = 0; s < mergeStart.size(); s++) {
                if (mergeStart.at(s).seq.size() < 20) {
                    continue;
                }

                if (s>=5) {
                    break;
                }

                for (int e = 0; e < mergeEnd.size(); e++) {


                    if (mergeEnd.at(e).seq.size() < 20) {
                        continue;
                    }

                    if (e>=5) {
                        break;
                    }

                    std::string seqSource = mergeStart.at(s).seq.substr(mergeStart.at(s).seq.size() - kmer);
                    std::string seqtarget = mergeEnd.at(e).seq.substr(e, kmer);

                    dynGraph.buildGraph(dynGraph.reverseString(mergeStart.at(s).seq));
                    dynGraph.buildGraph(dynGraph.reverseString(mergeEnd.at(e).seq));

                    for (std::string x:vUnMappedRead) {
                        dynGraph.buildGraph(dynGraph.reverseString(x));
                    }

                    GraphResult gr = dynGraph.getGraphResult(dynGraph.reverseString(seqSource),
                                                             dynGraph.reverseString(seqtarget));
                    gr = dynGraph.ReverseStringGraphResult(&gr);
                    if (50 < dynGraph.reverseString(gr.getMaxConcordant()).size()) {

//                        std::cout << "GR LEFT : " << gr.getMaxLeft() << std::endl;
//                        std::cout << "GR RIGTH : " << gr.getMaxRight() << std::endl;
//                        std::cout << "GR CONCORDANT : " << gr.getMaxConcordant() << std::endl;
//                        std::cout << "------------------" << std::endl;

                        BreakpointPosition tempBP;
                        tempBP.pos = n.getPosition();
                        tempBP.end = m.getPosition();
                        tempBP.frequency = mergeStart.at(s).count + mergeEnd.at(e).count;

                        tempBP.seqLeft = gr.getMaxLeft();
                        tempBP.seqRight = gr.getMaxRight();
                        tempBP.seqMatch = gr.getMaxConcordant();
                        tempBP.longmatch = gr.getMaxConcordant().size();
                        tempBP.score = n.getFrequency() + m.getFrequency();
                        tempBP.longmapstart = n.getLongMapping();
                        tempBP.longmapend = m.getLongMapping();

                        std::vector<uint8_t> tempmapq;
                        tempmapq.insert(tempmapq.end(), mergeStart.at(s).mapqlist.begin(),
                                        mergeStart.at(s).mapqlist.end());
                        tempmapq.insert(tempmapq.end(), mergeEnd.at(e).mapqlist.begin(), mergeEnd.at(e).mapqlist.end());

                        tempBP.mappingqualitylist = tempmapq;

//                        if (tempBP.frequency <= 2)
//                        {
//                            continue;
//                        }

                        added = true;
                        vectorBP.push_back(tempBP);

                    } else if (mergeStart.at(s).seq.size() < dynGraph.reverseString(gr.getMaxLeft()).size()) {

//                        std::cout << "GR LEFT : " << gr.getMaxLeft() << std::endl;
//                        std::cout << "GR RIGTH : " << gr.getMaxRight() << std::endl;
//                        std::cout << "GR CONCORDANT : " << gr.getMaxConcordant() << std::endl;
//                        std::cout << "------------------" << std::endl;

                        BreakpointPosition tempBP;
                        tempBP.pos = n.getPosition();
                        tempBP.end = m.getPosition();
                        tempBP.frequency = mergeStart.at(s).count + mergeEnd.at(e).count;

                        tempBP.seqLeft = gr.getMaxLeft();
                        tempBP.seqRight = gr.getMaxRight();
                        tempBP.seqMatch = gr.getMaxConcordant();
                        tempBP.longmatch = gr.getMaxConcordant().size();
                        tempBP.score = n.getFrequency() + m.getFrequency();
                        tempBP.longmapstart = n.getLongMapping();
                        tempBP.longmapend = m.getLongMapping();

                        std::vector<uint8_t> tempmapq;
                        tempmapq.insert(tempmapq.end(), mergeStart.at(s).mapqlist.begin(),
                                        mergeStart.at(s).mapqlist.end());
                        tempmapq.insert(tempmapq.end(), mergeEnd.at(e).mapqlist.begin(), mergeEnd.at(e).mapqlist.end());

                        tempBP.mappingqualitylist = tempmapq;

//                        if (tempBP.frequency <= 2)
//                        {
//                            continue;
//                        }

                        added = true;
                        vectorBP.push_back(tempBP);

                    } else if (mergeEnd.at(e).seq.size() < dynGraph.reverseString(gr.getMaxRight()).size()) {
//                        std::cout << "GR LEFT : " << gr.getMaxLeft() << std::endl;
//                        std::cout << "GR RIGTH : " << gr.getMaxRight() << std::endl;
//                        std::cout << "GR CONCORDANT : " << gr.getMaxConcordant() << std::endl;
//                        std::cout << "------------------" << std::endl;
//
                        BreakpointPosition tempBP;
                        tempBP.pos = n.getPosition();
                        tempBP.end = m.getPosition();
                        tempBP.frequency = mergeStart.at(s).count + mergeEnd.at(e).count;

                        tempBP.seqLeft = gr.getMaxLeft();
                        tempBP.seqRight = gr.getMaxRight();
                        tempBP.seqMatch = gr.getMaxConcordant();
                        tempBP.longmatch = gr.getMaxConcordant().size();

                        tempBP.score = n.getFrequency() + m.getFrequency();
                        tempBP.longmapstart = n.getLongMapping();
                        tempBP.longmapend = m.getLongMapping();

                        std::vector<uint8_t> tempmapq;
                        tempmapq.insert(tempmapq.end(), mergeStart.at(s).mapqlist.begin(),
                                        mergeStart.at(s).mapqlist.end());
                        tempmapq.insert(tempmapq.end(), mergeEnd.at(e).mapqlist.begin(), mergeEnd.at(e).mapqlist.end());

                        tempBP.mappingqualitylist = tempmapq;

//                        if (tempBP.frequency <= 2)
//                        {
//                            continue;
//                        }

                        added = true;
                        vectorBP.push_back(tempBP);
                    }
                }
            }

        }

    }

}

bool RefineInsertion::checkBetween(int32_t pos, int32_t targetPos, int32_t minusoverlapped, int32_t plusoverlapped) {
    if (targetPos + minusoverlapped > pos) {
        return false;
    }

    if (targetPos + plusoverlapped < pos) {
        return false;
    }

    return true;
}

void RefineInsertion::convertMapSC() {
    vectorSCStart = convertMapSCToVector(mapSCStart);
    std::sort(vectorSCStart.begin(), vectorSCStart.end());

    vectorSCEnd = convertMapSCToVector(mapSCEnd);
    std::sort(vectorSCEnd.begin(), vectorSCEnd.end());
}

std::vector<InsertionPositionDetail>
RefineInsertion::convertMapSCToVector(std::map<int32_t, InsertionPositionDetail> mapSC) {

    std::vector<InsertionPositionDetail> vectorSC;
    int count = 0;
    for (auto x : mapSC) {

        if (x.second.getLongMapping() < 30) {
            continue;
        }

        vectorSC.push_back(x.second);

        count++;
    }

    return vectorSC;
}

void RefineInsertion::clearMapSC() {
    mapSCStart.clear();
    mapSCEnd.clear();
}

std::vector<RefineInsertion::CountRefineSeq>
RefineInsertion::mergeString(InsertionPositionDetail fragmentlist, bool fromstart) {

    std::vector<CountRefineSeq> tempSeq;

    int count = 0;
    for (std::string n : fragmentlist.getSeqList())
        //  for (int in=0;in < fragmentlist.)
    {

        bool added = false;
        for (int i = 0; i < tempSeq.size(); i++) {
            if (compareEditDistance(n, tempSeq.at(i).seq, fromstart)) {
                if (n.size() > tempSeq.at(i).seq.size()) {
                    tempSeq.at(i).seq = n;
                    tempSeq.at(i).count++;
                    tempSeq.at(i).mapqlist.push_back(fragmentlist.getMapQList().at(count));
                }
                added = true;
                // break;
            }
        }

        if (!added) {
            CountRefineSeq tempCRS;
            tempCRS.seq = n;
            tempCRS.count++;
            tempCRS.mapqlist.push_back(fragmentlist.getMapQList().at(count));
            tempSeq.push_back(tempCRS);
        }
        count++;
    }

    return tempSeq;
}

bool RefineInsertion::compareEditDistance(std::string s1, std::string s2, bool fromstart) {
    std::string temps1 = s1;
    std::string temps2 = s2;
    substringSeq(&temps1, &temps2, fromstart);

    EditDistance editdistance;
    int editpoint = editdistance.Compare(&temps1, &temps2);

    if (editpoint < 3) {
        return true;
    }

    return false;
}

void RefineInsertion::substringSeq(std::string *s1, std::string *s2, bool fromstart) {
    std::string temps1;
    std::string temps2;
    int diffsize = s1->size() - s2->size();
    if (s1->size() > s2->size()) {
        if (!fromstart) {
            temps1 = s1->substr(s1->length() - s2->size());
            temps2 = *s2;
        } else {
            temps1 = s1->substr(0, s2->size());
            temps2 = *s2;
        }
    } else {
        if (!fromstart) {
            temps2 = s2->substr(s2->length() - s1->size());
            temps1 = *s1;
        } else {
            temps2 = s2->substr(0, s1->size());
            temps1 = *s1;
        }
    }

    *s1 = temps1;
    *s2 = temps2;
}

bool RefineInsertion::getOverlappedSeq(std::vector<CountRefineSeq> startSeq, std::vector<CountRefineSeq> endSeq,
                                       int *frequency, int *longmatch, std::vector<uint8_t> *mapq, std::string *seq1,
                                       std::string *seq2) {

    for (CountRefineSeq n : startSeq) {
        if (n.seq.size() < 30) {
            continue;
        }


        for (CountRefineSeq m : endSeq) {
            if (m.seq.size() < 30) {
                continue;
            }

            if (n.count + m.count <= 2) {
                continue;
            }

            if ((evidence.getMark() == "")) {
                SmithWaterman swm(&n.seq, 0, false);
                int maxmatch = swm.findMaxMatchInsertion(&m.seq);

                int seq1MatchSize = n.seq.size() - maxmatch;
                int seq2MatchSize = m.seq.size() - maxmatch;

                if (seq1MatchSize + seq2MatchSize + maxmatch < 50) {
                    continue;
                }

                if (maxmatch <= getDivider(samplestat->getReadLength(), 15, 100, 15)) {
                    continue;
                }

                *frequency = n.count + m.count;
                std::vector<uint8_t> tempmapq;
                tempmapq.insert(tempmapq.end(), n.mapqlist.begin(), n.mapqlist.end());
                tempmapq.insert(tempmapq.end(), m.mapqlist.begin(), m.mapqlist.end());
                *mapq = tempmapq;
                *longmatch = maxmatch;

                // if (samplestat->get)

                if (maxmatch > 15) {
                    return true;
                }
            } else {

                if (n.seq.size() <= getDivider(samplestat->getReadLength(), 20, 100, 25)) {
                    continue;
                }

                if (m.seq.size() <= getDivider(samplestat->getReadLength(), 20, 100, 25)) {
                    continue;
                }

                *frequency = n.count + m.count;
                std::vector<uint8_t> tempmapq;
                tempmapq.insert(tempmapq.end(), n.mapqlist.begin(), n.mapqlist.end());
                tempmapq.insert(tempmapq.end(), m.mapqlist.begin(), m.mapqlist.end());
                *mapq = tempmapq;
                if (m.seq.size() >= n.seq.size()) {
                    *longmatch = m.seq.size();
                } else {
                    *longmatch = n.seq.size();
                }

                return true;
            }
        }
    }

    return false;
}
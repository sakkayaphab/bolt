#include "specifyevidence.h"

SpecifyEvidence::SpecifyEvidence()
{
    setMaxRemainBufferWriteEvidenceFile(10);
}

void SpecifyEvidence::setRead(bam1_t *bam_read, bam_hdr_t *m_bam_header)
{
    read = bam_read;
    bam_header = m_bam_header;

    readparser.setBamRead(read);
    readparser.setBamHeader(bam_header);
}

void SpecifyEvidence::setSampleStat(SampleStat *m_samplestat)
{
    samplestat = m_samplestat;
}

void SpecifyEvidence::showSizeFinalEvidence()
{
    std::cout << "--------+ " << svtype << " +--------" << std::endl;
    std::cout << "size evidence : " << finalEvidence.size() << std::endl;
}

bool SpecifyEvidence::checkBetween(int32_t pos, int32_t targetPos, int32_t overlapped)
{
    if (targetPos - overlapped > pos)
    {
        return false;
    }

    if (targetPos + overlapped < pos)
    {
        return false;
    }

    return true;
}

void SpecifyEvidence::setReadDepthHelper(ReadDepthHelper *m_readdepthHelper) {
    readdepthHelper = m_readdepthHelper;
}

int SpecifyEvidence::findOverlappedOnlyPos(uint32_t overlapped, uint32_t pos)
{
    for (int i = 0; i < preCollectSV.size(); i++)
    {
        if (checkBetween(pos, preCollectSV.at(i).getPosDiscordantRead(), overlapped))
        {
            return i;
        }

        if (pos > preCollectSV.at(i).getPosDiscordantRead() + overlapped)
        {
            return -1;
        }
    }
    return -1;
}

int SpecifyEvidence::findOverlapped(int32_t overlapped, int32_t pos, int32_t mpos)
{
    for (int i = 0; i < preCollectSV.size(); i++)
    {
        if (checkBetween(pos, preCollectSV.at(i).getPosDiscordantRead(), overlapped) && checkBetween(mpos,
                                                                                                     preCollectSV.at(i).getEndDiscordantRead(),
                                                                                                     overlapped))
        {
            return i;
        }

        if (pos > preCollectSV.at(i).getPosDiscordantRead() + overlapped)
        {
            return -1;
        }
    }
    return -1;
}

void SpecifyEvidence::showAllFinalEvidence()
{
    for (auto n : finalEvidence)
    {
        int idiff = n.getEndDiscordantRead() - n.getLastPosDiscordantRead();
        // uint32_t diff = std::abs(idiff);
        std::cout << "Pos :"
                  << n.getPosDiscordantRead()
                  << ", Last Pos :"
                  << n.getLastPosDiscordantRead()
                  << ", Mate Pos :"
                  << n.getEndDiscordantRead()
                  << ", Last Mate Pos :"
                  << n.getLastEndDiscordantRead()
                  << ", Frequency :"
                  << n.getFrequency()
                  << ", Forward Direction :"
                  << n.getForwardDirection()
                  << ", Diff Pos :"
                  << idiff
                  << std::endl;
    }
}

void SpecifyEvidence::setOutputPath(std::string path)
{
    outputpath = path;
}

void SpecifyEvidence::writeFinalEvidenceAndClear()
{
    std::ofstream myfile;
    myfile.open(outputpath, std::ios_base::app);
    if (myfile.is_open())
    {
        for (auto n : finalEvidence)
        {
            myfile << n.convertToVcfString()
                   << std::endl;
        }
        finalEvidence.clear();

        myfile.close();
    }
}

int SpecifyEvidence::getMaxRemainBufferWriteEvidenceFile() const {
    return maxRemainBufferWriteEvidenceFile;
}

void SpecifyEvidence::setMaxRemainBufferWriteEvidenceFile(int32_t bufferWriteEvidenceFile) {
    SpecifyEvidence::maxRemainBufferWriteEvidenceFile = bufferWriteEvidenceFile;
}

void SpecifyEvidence::writeBufferEvidenceFile()
{
    if (finalEvidence.size() > getMaxRemainBufferWriteEvidenceFile())
    {
        writeFinalEvidenceAndClear();
    }
}
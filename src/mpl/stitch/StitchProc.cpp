/*************************************************************************
> File Name: StitchProc.h
> Author: Qi Sun
> Mail: qsun@cse.cuhk.edu.hk
> Created Time: 27/08/2018
************************************************************************/

#include "StitchProc.h"

SIMPLEMPL_BEGIN_NAMESPACE

StitchProc::StitchProc()
{
	this->reset(true);
}

StitchProc::StitchProc(double pitch)
{
	this->reset(true);
	this->PITCH = pitch;
}

StitchProc::~StitchProc()
{
	if (m_db) delete m_db;
}


void StitchProc::read_cmd(int argc, char** argv)
{
    ControlParameter tmpParms;
    CmdParser cmd(tmpParms);
    mplAssertMsg(cmd(argc,argv), "failed to parse command");
    
    if (tmpParms.shape_mode==ShapeModeEnum::RECTANGLE)
        m_db = new LayoutDBRect;
    else
        m_db = new LayoutDBPolygon;
    tmpParms.swap(m_db->parms);
}

void StitchProc::reset(bool init)
{
	// release memory and set to initial value
	if (!init)
	{
		if (m_db) delete m_db;
		std::vector<uint32_t>().swap(m_vVertexOrder);
		std::vector<std::vector<uint32_t>>().swap(m_Touch);
		std::vector<std::vector<uint32_t>>().swap(m_DPL);
		std::vector<uint32_t>().swap(m_vCompId);
		std::vector<std::pair<uint32_t, uint32_t>>().swap(m_vConflict);
	}
	m_db = NULL;
	m_comp_cnt = 0;
}

void StitchProc::read_gds()
{
	char buf[256];
	mplSPrint(kINFO, buf, "reading input files takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer(buf);
	mplPrint(kINFO, "Reading input file %s\n", m_db->input_gds().c_str());
	// read input gds file 
	GdsReader reader(*m_db);
	mplAssertMsg(reader(m_db->input_gds()), "failed to read %s", m_db->input_gds().c_str());
	// must call initialize after reading to initize the rtree.
	m_db->initialize_data();

	// report data 
	m_db->report_data();
}

void StitchProc::write_gds()
{
	char buf[256];
	mplSPrint(kINFO, buf, "writing output file takes %%t seconds CPU, %%w seconds real\n");
	boost::timer::auto_cpu_timer timer(buf);
	if (m_db->output_gds().empty())
	{
		mplPrint(kWARN, "Output file not specified, no file generated\n");
		return;
	}
	// write output gds file 
	GdsWriter writer;
	mplPrint(kINFO, "Write output gds file: %s\n", m_db->output_gds().c_str());
	writer(m_db->output_gds(), *m_db, m_vConflict, m_DPL, m_db->strname, m_db->unit*1e+6);
}

/*
void StitchProc::solve()
{
	char buf[256];
	boost::timer::auto_cpu_timer timer(buf);

	if (m_db->vPatternBbox.empty())
	{
		mplPrint(kWARN, "No patterns found in specified layers.\n");
		return;
	}
	this->construct_graph();

}
*/

// This method is only for rectangles.
void StitchProc::relation4rect()
{
    uint32_t vertex_num = m_db->vPatternBbox.size();
    m_DPL.resize(vertex_num);
    for (int i =0; i<m_db->vPatternBbox.size();i++)
    {
        rectangle_pointer_type const& pPattern = m_db->vPatternBbox[i];

		m_layout_left = (m_layout_left > gtl::xl(*pPattern)) ? gtl::xl(*pPattern) : m_layout_left;
		m_layout_right = (m_layout_right < gtl::xh(*pPattern)) ? gtl::xh(*pPattern) : m_layout_right;
		m_layout_top = (m_layout_top < gtl::yh(*pPattern)) ? gtl::yh(*pPattern) : m_layout_top;
		m_layout_down = (m_layout_down > gtl::yl(*pPattern)) ? gtl::yl(*pPattern) : m_layout_down;


        std::vector<uint32_t>& vNeiDPL = m_DPL[i];
		std::vector<uint32_t>& vNeiSafe = m_Safe[i];
        
        rectangle_type rect(*pPattern);
		gtl::bloat(rect, gtl::HORIZONTAL, m_db->coloring_distance);
		gtl::bloat(rect, gtl::VERTICAL, m_db->coloring_distance);
        for (rtree_type::const_query_iterator itq=m_db->tPatternBbox.qbegin(bgi::intersects(rect));
             itq!=m_db->tPatternBbox.qend(); ++itq) 
		{
            rectangle_pointer_type const& pAdjPattern = *itq;
			if (pAdjPattern != pPattern)
			{
				mplAssert(pAdjPattern->pattern_id() != pPattern->pattern_id());
				coordinate_difference distance = m_db->euclidean_distance(*pAdjPattern, *pPattern);
				if (distance < m_db->coloring_distance)
					vNeiDPL.push_back(pAdjPattern->pattern_id());
				else if (distance < m_db->coloring_distance + PITCH)
					vNeiSafe.push_back(pAdjPattern->pattern_id());
			}
        }



		vNeiDPL.swap(vNeiDPL);
		vNeiSafe.swap(vNeiSafe);
    }
    return;
}


void StitchProc::projection()
{
	
	uint32_t vertex_num = m_db->vPatternBbox.size();
	m_interRect.resize(vertex_num);
	// traverse all vertices in the layout graph
	for (int i = 0; i < vertex_num; i++)
	{
		// step 1 : init vinterRect, store all overlapping bounding boxes
		rectangle_pointer_type const& pPattern = m_db->vPatternBbox[i];
		bool isHor = whetherHor(pPattern);
		std::vector<rectangle_pointer_type> vinterRect;
		vinterRect.clear();
		std::vector<uint32_t>::iterator itr = m_DPL[pPattern->pattern_id()].begin();
		std::vector<uint32_t>::iterator end = m_DPL[pPattern->pattern_id()].end();
		for (; itr != end; itr++)
		{
			rectangle_pointer_type AdjPattern = m_db->vPatternBbox[*itr];
			rectangle_type extend (gtl::xl(*AdjPattern) - MINDPLDIST, gtl::yl(*AdjPattern) - MINDPLDIST, gtl::xh(*AdjPattern) - MINDPLDIST, gtl::yh(*AdjPattern) - MINDPLDIST);
			std::deque<rectangle_type> output;
			bg::intersection(extend, pPattern, output);
			BOOST_FOREACH(rectangle_type p, output)
			{
				vinterRect.push_back(&p);
			}
		} // for itr, present Pattern's neighbors

		// step 2 : init vstitches, generate all the candidate stitches for present Pattern
		std::vector <coordinate_type> vstitches;
		rectangle_pointer_type tempRect = new rectangle_type(gtl::xl(*pPattern), gtl::yl(*pPattern), gtl::xh(*pPattern), gtl::yh(*pPattern));
		coordinate_type left = 0, right = 0;
		if (isHor) { left = gtl::xl(*tempRect); right = gtl::xh(*tempRect); }
		else { left = gtl::yl(*tempRect); right = gtl::yh(*tempRect); }

		if (true == dplstitch)
			StitchGenerateDPL_Points(tempRect, vinterRect, vstitches, left, right, false);
		else
			StitchGenerateTPL_Points(tempRect, vinterRect, vstitches, left, right, false);

		// Step 3 : clean up vinterRect and tempRect, vstitches should be preserved
		// delete all the temperal rectangles
		for (int j = 0; j < vinterRect.size(); j++)
			delete vinterRect[j];
		delete tempRect;

		// Step 4 : check the positions' legalities
		if(isHor)
		{
			for(int j=0;j<vstitches.size();j++)
			{
				int pos = vstitches[i];
				assert(pos>gtl::xl(*pPattern) && pos < gtl::xh(*pPattern));
			} // for j
		}
		else
		{
			for(int j=0; j< vstitches.size();j++)
			{
				int pos = vstitches[j];
				assert(pos>gtl::yl(*pPattern) && pos < gtl::yh(*pPattern));
			} // for j
		} // for if

		// Step 5 : generate the final stitch candidates
		if(vstitches.size() <= 0)	// have not been split
		{
			m_final_rectangle.push_back(*pPattern);
		}
		else 
		{
			int size = vstitches.size();
			for(int pos=0; pos < vstitches.size(); pos ++)
			{
				rectangle_pointer_type newRect = new rectangle_type;
				if (pos == 0)
				{
					if (isHor) {
						newRect->
					}
				}
			}
		}
	}
}


bool StitchProc::whetherHor(rectangle_pointer_type tmp)
{
	double xl = gtl::xl(*tmp);
	double yl = gtl::yl(*tmp);
	double xh = gtl::xh(*tmp);
	double yh = gtl::yh(*tmp);
	return (xh - xl) > (yh - yl);
}


void StitchProc::StitchGenerateDPL_Points(const rectangle_pointer_type pRect, const std::vector<rectangle_pointer_type> vinterRect, std::vector <coordinate_type> vstitches, coordinate_type left, coordinate_type right, bool vddstitch = false)
{
	vstitches.clear();
	bool isHor = whetherHor(pRect);

	// No stitch in VDD/GND
	if (false == vddstitch)
	{
		if (isHor && getWidth(pRect) > (m_layout_right - m_layout_left)*0.6) return;
		if (!isHor && getHeight(pRect) > (m_layout_top - m_layout_down)*0.6) return;
	}

	// tempset stores all the possible locations.
	std::set<int> tempset;
	tempset.insert(left);
	tempset.insert(right);

	// travers all the conflicts
	for (int i = 0; i < vinterRect.size(); i++)
	{
		if (isHor)
		{
			tempset.insert(gtl::xl(*vinterRect[i]));
			tempset.insert(gtl::xh(*vinterRect[i]));
		}
		else
		{
			tempset.insert(gtl::yl(*vinterRect[i]));
			tempset.insert(gtl::yh(*vinterRect[i]));
		}
	} // for i 

	// tempvec is used to store the possible stitches in order.
	std::vector<int> tempvec;
	std::set<int>::iterator itr = tempset.begin();
	std::set<int>::iterator end = tempset.end();

	for (; itr != end; itr++) tempvec.push_back(*itr);
	std::sort(tempvec.begin(), tempvec.end());

	// accumulate the overlapping interseciton
	std::vector<std::pair<std::pair<int, int>, int>> vstages;
	for (int i = 1; i < tempvec.size(); i++)
		vstages.push_back(std::make_pair(std::make_pair(tempvec[i - 1], tempvec[i]), 0));
	for (int i = 0; i < vstages.size(); i++)
	{
		int lpos = vstages[i].first.first;
		int rpos = vstages[i].first.second;
		int cumu = 0;
		for (int j = 0; j < vinterRect.size(); j++)
		{
			if (isHor) {
				if (lpos < gtl::xl(*vinterRect[j])) continue;
				if (rpos > gtl::xh(*vinterRect[j])) continue;
				cumu++;
			}
			else
			{
				if (lpos < gtl::yl(*vinterRect[j])) continue;
				if (rpos > gtl::yh(*vinterRect[j])) continue;
				cumu++;
			}
		} // for j
		vstages[i].second = cumu;
	} // for i

	// find zeros in projection sequence
	for (int i = 1; i < vstages.size() - 1; i++)
	{
		if (vstages[i].second > 0) continue;
		int pos = (vstages[i].first.first + vstages[i].first.second) / 2;
		vstitches.push_back(pos);
	}
	std::sort(vstitches.begin(), vstitches.end());
}


void StitchProc::StitchGenerateTPL_Points(const rectangle_pointer_type pRect, const std::vector<rectangle_pointer_type> vinterRect, std::vector <coordinate_type> vstitches, coordinate_type left, coordinate_type right, bool vddstitch = false)
{
	vstitches.clear();
	bool isHor = whetherHor(pRect);

	// use tempset to store all possible stitches
	// then sort the stitches and store them into tempvec
	std::set<int> tempset;
	// traverse all the conflicts to generate all possible stitches
	tempset.insert(left);
	tempset.insert(right);
	for (int i = 0; i < vinterRect.size(); i++)
	{
		if (isHor)
		{
			tempset.insert(gtl::xl(*vinterRect[i]));
			tempset.insert(gtl::xh(*vinterRect[i]));
		}
		else
		{
			tempset.insert(gtl::yl(*vinterRect[i]));
			tempset.insert(gtl::yh(*vinterRect[i]));
		}
	} // for i
	
	std::vector<int> tempvec;
	std::set<int>::iterator itr = tempset.begin(), end = tempset.end();
	for (; itr != end; itr++)	tempvec.push_back(*itr);
	std::sort(tempvec.begin(), tempvec.end());

	// Step 1 : init vstages
	// all the points are stored at tmpvec
	std::vector<std::pair<std::pair<int, int>, int>> vstages;
	for (int i = 1; i < tempvec.size(); i++)
		vstages.push_back(std::make_pair(std::make_pair(tempvec[i - 1], tempvec[i]), 0));
	for (int i = 0; i < vstages.size(); i++)
	{
		int lpos = vstages[i].first.first;
		int rpos = vstages[i].first.second;
		int cumu = 0;
		for (int j = 0; j < vinterRect.size(); j++)
		{
			if (isHor) {
				if (lpos < gtl::xl(*vinterRect[j])) continue;
				if (rpos > gtl::xh(*vinterRect[j])) continue;
				cumu++;
			}
			else
			{
				if (lpos < gtl::yl(*vinterRect[j])) continue;
				if (rpos > gtl::yh(*vinterRect[j])) continue;
				cumu++;
			}
		} // for j
		vstages[i].second = cumu;
	} // for i

	// Step 2 : 
	// copy into projection sequence and add default terminal zeros
	std::vector<std::pair<std::pair<int, int>, int>> vSequence;
	if (vstages.size() <= 0) return;
	if (vstages[0].second != 0)
	{
		vSequence.push_back(std::make_pair(std::make_pair(left, left), 0));
	}
	int size = vstages.size();
	for (int i = 0; i < size; i++)
	{
		vSequence.push_back(vstages[i]);
	}
	if (vstages[size - 1].second != 0)
	{
		vSequence.push_back(std::make_pair(std::make_pair(right, right), 0));
	}

	// Step 3 : 
	// find zeros in vSequence[].second
	std::vector<int> vzeroids;
	for (int i = 0; i < vSequence.size(); i++)
	{
		if (vSequence[i].second > 0) continue;
		vzeroids.push_back(i);
	}

	// Step 4 : 
	// all return stitch positions are in vstitches
	// These operations are a little confusing.
	for (int i = 0; i < vzeroids.size() - 1; i++)
	{
		int pos1 = vzeroids[i];
		int pos2 = vzeroids[i + 1];
		if (0 == i) assert(0 == pos1);
		else if (2 == pos1)		// remove the useless stitch
		{
			bool bFind = true;
			if (vSequence.size() < 5) bFind = false;
			else if (i != 1) bFind = false;
			else if (1 != vSequence[1].second) bFind = false;
			else if (1 != vSequence[3].second) bFind = false;
			else if (0 != vSequence[4].second) bFind = false;
			if (bFind == true) continue;
			int pos = (vSequence[2].first.first + vSequence[2].first.second) / 2;
			vstitches.push_back(pos);
		}
		else if (pos1 == vSequence.size() - 3)
		{
			assert(i == vzeroids.size() - 2);
			bool bFind = true;
			int zsize = vzeroids.size();
			if (vSequence.size() < 5) bFind = false;
			else if (i != zsize - 2) bFind = false;
			else if (1 != vSequence[pos1 + 1].second) bFind = false;
			else if (1 != vSequence[pos1 - 1].second) bFind = false;
			else if (0 != vSequence[pos1 - 2].second) bFind = false;
			if (bFind == true) continue;
			int pos = (vSequence[pos1].first.first + vSequence[pos1].first.second) / 2;
			vstitches.push_back(pos);
		}
		else
		{
			int pos = (vSequence[pos1].first.first + vSequence[pos1].first.second) / 2;
			vstitches.push_back(pos);
		}

		// search lost stitch in vSequence[pos1 --> pos2]
		double maxValue = 0.9;
		int posLost = pos1 + 2;
		if (pos2 - pos1 < 4) continue;
		for (int i = pos1 + 2; i < pos2 - 1; i++)
		{
			if (vSequence[i - 1].second <= vSequence[i].second) continue;
			if (vSequence[i + 1].second <= vSequence[i].second) continue;
			int mind = std::min(vSequence[i - 1].second - vSequence[i].second, vSequence[i + 1].second - vSequence[i].second);
			int diff = std::abs(vSequence[i + 1].second - vSequence[i - 1].second);
			double value = (double)mind + (double)diff * 0.1;
			if (value > maxValue)
			{
				maxValue = value;
				posLost = i;
			}
		}
		if (maxValue > 0.9)
		{
			int pos = (vSequence[posLost].first.first + vSequence[posLost].first.second) / 2;
			vstitches.push_back(pos);
		}
	}
	std::sort(vstitches.begin(), vstitches.end());
}



SIMPLEMPL_END_NAMESPACE
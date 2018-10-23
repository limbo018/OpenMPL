void SimpleMPL::GenerateStitchPositionBei(const rectangle_pointer_type pRect,
	const std::vector<rectangle_type> vinterRect,
	std::vector <coordinate_type> vstitches, const coordinate_type lower,
	const coordinate_type upper)
{
#ifdef QDEBUG
	std::cout << "pattern id : " << pRect->pattern_id() << "\tlower : " << lower << "\tupper : " << upper << std::endl;
#endif
	// ================================================================================
	// step 1 : generate candidate stitches' positions according to the intersections
	// ================================================================================
	bool isHor = whetherHorizontal(pRect);
	std::set<coordinate_type> tempSet;
	tempSet.insert(lower);
	tempSet.insert(upper);
	for (std::vector<rectangle_type>::const_iterator it = vinterRect.begin(); it != vinterRect.end(); it++)
	{
		if (isHor)
		{
			tempSet.insert(gtl::xl(*it));
			tempSet.insert(gtl::xh(*it));
		}
		else
		{
			tempSet.insert(gtl::yl(*it));
			tempSet.insert(gtl::yh(*it));
		}
	}
	std::vector<coordinate_type> tempVec;
	for(std::set<coordinate_type>::iterator itr = tempSet.begin(); itr != tempSet.end(); itr++)
		tempVec.push_back(*itr);
	// sort all the positions
	sort(tempVec.begin(), tempVec.end());
#ifdef QDEBUG
	// ouput tempSet
	std::cout << "==== tempSet ====" << std::endl; 
	for(std::set<coordinate_type>::iterator itr = tempSet.begin(); itr != tempSet.end(); itr++)
		std::cout << *itr << std::endl;
	// output tempvec
	std::cout << "==== tempVec ====" << std::endl;
	for(uint32_t i = 0; i < tempVec.size(); i ++)
	{
		std::cout << i << " : \t" << tempVec[i] << std::endl; 
	}
#endif
	// ================================================================================
	// step 2 : generate stages according to the stitch positions
	//			that is each pairwise of neighboring positions generates a stage.
	// ================================================================================
	std::vector<std::pair<std::pair<coordinate_type, coordinate_type>, uint32_t> > vStages;
	for (uint32_t i = 1; i < tempVec.size(); i++)
		vStages.push_back(std::make_pair(std::make_pair(tempVec[i - 1], tempVec[i]), 0));
	// calculate the times every stage covered by all the intersections.
	for (std::vector<std::pair<std::pair<coordinate_type, coordinate_type>, uint32_t> >::iterator it = vStages.begin();
		it != vStages.end(); it++)
	{
		coordinate_type left = it->first.first;
		coordinate_type right = it->first.second;
		int overlapping_count = 0;
		for (std::vector<rectangle_type>::const_iterator itInt = vinterRect.begin(); itInt != vinterRect.end(); itInt++)
		{
			if (isHor)
			{
				if (left < gtl::xl(*itInt)) continue;
				if (right > gtl::xh(*itInt)) continue;
				overlapping_count++;
			}
			else
			{
				if (left < gtl::yl(*itInt)) continue;
				if (right > gtl::yh(*itInt)) continue;
				overlapping_count++;
			}
		}
		it->second = overlapping_count;
	}

	// ================================================================================
	// step 3 : add default terminal zeros
	//			This will be used in the next step.
	// ================================================================================
	if (vStages.size() <= 0) return;
	if (vStages[0].second != 0)
		vStages.insert(vStages.begin(), std::make_pair(std::make_pair(lower, lower), 0));
	if (vStages[vStages.size() - 1].second != 0)
		vStages.push_back(std::make_pair(std::make_pair(upper, upper), 0));
#ifdef DEBUG
	
	std::cout << "DEBUG_PROJECTION| vStages = ";
	for (std::vector<std::pair<std::pair<coordinate_type, coordinate_type>, uint32_t> >::iterator it = vStages.begin();
		it != vStages.end(); it++)
		std::cout << it->second;
	std::cout << std::endl;
	
#endif

	// ================================================================================
	// step 4: find the stages with zero overlapping_count
	//		   The stitches will be chosen from these stages.
	// ================================================================================
	std::vector<uint32_t> vZeroIds;
	for (uint32_t i = 0; i < vStages.size(); i++)
	{
#ifdef QDEBUG
		std::cout << "vStages[" << i <<  "] \t" << vStages[i].first.first << " -- " << vStages[i].first.second << "\t count = " << vStages[i].second << std::endl;
#endif
		if (vStages[i].second > 0) continue;
		vZeroIds.push_back(i);
	}
#ifdef QDEBUG
	// output vZeroIds
	std::cout << "Zero stages : " << std::endl;
	for (uint32_t i = 0; i < vZeroIds.size(); i++ )
	{
		std::cout << i << " : " << "vStages[" << vZeroIds[i] << "] " << vStages[vZeroIds[i]].first.first << " -- " << vStages[vZeroIds[i]].first.second << std::endl;
	}
	std::cout << "==== Choose stitches ==== " << std::endl << "vStages.size() :  " << vStages.size() << std::endl;
#endif
	// ================================================================================
	// step 5: choose stitches from vZeroIds
	// ================================================================================
	std::vector<coordinate_type>().swap(vstitches);
	// The operations here are very confusing.
	for (uint32_t i = 0; i < vZeroIds.size() - 1; i++)
	{
		uint32_t pos1 = vZeroIds[i];
		uint32_t pos2 = vZeroIds[i + 1];
#ifdef QDEBUG
		std::cout << "i = " << i << " \t" << "pos1 = " << pos1 << " \t" << "pos2 = " << pos2 << std::endl;
#endif
		// since ((lower, lower), 0) has been added into vStages, so pos1 must be 0.
		if (i == 0) mplAssertMsg(0 == pos1, "pos1 %d doesn't equal to 0", pos1);
		// remove the useless stitches
		else if (pos1 == 2)
		{
			bool find = false;
			if (vStages.size() < 5) find = true;
			else if (i != 1) find = true;
			else if (1 != vStages[1].second) find = true;
			else if (1 != vStages[3].second) find = true;
			else if (0 != vStages[4].second) find = true;
			if (find == false)
				continue;
			coordinate_type position = (vStages[2].first.first + vStages[2].first.second) / 2;
			vstitches.push_back(position);
		}
		else if (pos1 == vStages.size() - 3)
		{
			std::cout << "i = " << i << "\t" << " vZeroIds.size() = " << vZeroIds.size() << std::endl;
			mplAssert(i == vZeroIds.size() - 2);
			bool find = false;
			uint32_t zsize = vZeroIds.size();
			if (vStages.size() < 5) find = true;
			else if (i != zsize - 2) find = true;
			else if (1 != vStages[pos1 + 1].second) find = true;
			else if (1 != vStages[pos1 - 1].second) find = true;
			else if (0 != vStages[pos1 - 2].second) find = true;
			if (find == false) continue;
			coordinate_type position = (vStages[pos1].first.first + vStages[pos1].first.second) / 2;
			vstitches.push_back(position);
		}
		else
		{
			coordinate_type position = (vStages[pos1].first.first + vStages[pos1].first.second) / 2;
			vstitches.push_back(position);
		}

		// search lost stitch in vStages[pos1 --> pos2]
		double maxValue = 0.9;
		uint32_t posLost = pos1 + 2;
		if (pos2 - pos1 < 4) continue;
		for (uint32_t i = pos1 + 2; i < pos2 - 1; i++)
		{
			if (vStages[i - 1].second <= vStages[i].second) continue;
			if (vStages[i + 1].second <= vStages[i].second) continue;
			uint32_t mind = std::min(vStages[i - 1].second - vStages[i].second, vStages[i + 1].second - vStages[i].second);
			uint32_t diff = std::abs(static_cast<int>(vStages[i + 1].second - vStages[i - 1].second));
			double value = (double)mind + (double)diff * 0.1;
			if (value > maxValue)
			{
				maxValue = value;
				posLost = i;
			}
		}
		if (maxValue > 0.9)
		{
			uint32_t position = (vStages[posLost].first.first + vStages[posLost].first.second) / 2;
			vstitches.push_back(position);
		}
	}
	sort(vstitches.begin(), vstitches.end());
#ifdef DEBUG
	
	 std::cout << "DEBUG| output the vStages: " << std::endl;
	for (std::vector<std::pair<std::pair<coordinate_type, coordinate_type>, uint32_t> >::iterator it = vStages.begin();
		it != vStages.end(); it++)
		std::cout << it->second << "[" << it->first.first << ", " << it->first.second << "]" << std::endl;
	std::cout << std::endl;
	
#endif

#ifdef QDEBUG
	std::cout << "DEBUG| output the vstitches: " << std::endl;
	std::cout << "lower = " << lower << std::endl;
	std::cout << "upper = " << upper << std::endl;
	std::cout << 
	for (std::vector<coordinate_type>::iterator it = vstitches.begin(); it != vstitches.end(); it++)
		std::cout << *it << " ";
	std::cout << std::endl;
#endif
	return;
}
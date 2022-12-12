#pragma once

#define forAll(tempvector) \
    for(unsigned int i=0; i<tempvector.size(); i++) \
        for(unsigned int j=0; j<tempvector[i].size(); j++)
#define forAllInternal(tempvector) \
    for(unsigned int i=1; i<tempvector.size()-1; i++) \
        for(unsigned int j=1; j<tempvector[i].size()-1; j++)
#define forAllInternalUCVs(tempvector) \
       for(unsigned int i=1; i<tempvector.size()-2; i++) \
        for(unsigned int j=1; j<tempvector[i].size()-1; j++)
#define forAllInternalVCVs(tempvector) \
       for(unsigned int j=1; j<tempvector.size()-2; j++) \
        for(unsigned int i=1; i<tempvector[j].size()-1; i++)
#define forAllInternalX(tempvector) \
    for(unsigned int i=1; i<tempvector.size()-1; i++) \
        for(unsigned int j=1; j<tempvector[i].size(); j++)
#define forAllInternalY(tempvector) \
    for(unsigned int i=1; i<tempvector.size(); i++) \
        for(unsigned int j=1; j<tempvector[i].size()-1; j++)
#define forSouthBoundary(tempvector) \
    for(unsigned int i=1; i<tempvector.size()-1; i++) \
        for(unsigned int j=0; j<1; j++)
#define forNorthBoundary(tempvector) \
    for(unsigned int i=1; i<tempvector.size()-1; i++) \
        for(unsigned int j=tempvector[i].size()-1; j<tempvector[i].size(); j++)
#define forWestBoundary(tempvector) \
  for(unsigned int j=1; j<tempvector[0].size(); j++) \
    for(unsigned int i=1; i<2; i++)
#define forEastBoundary(tempvector) \
  for(unsigned int j=1; j<tempvector[0].size(); j++) \
        for(unsigned int i=tempvector[j].size()-2; i<tempvector[j].size()-1; i++)


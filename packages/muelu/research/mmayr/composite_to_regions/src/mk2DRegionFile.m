function [globalDims,localDims,relCorner,regionCorner] = mk2RegionFile(filename)
%
%  writes mapping to a file based on information below.
%

if (strcmp(filename, 'caseFive') == true)  % caseFive
% case 5 : regions may cross processors but no processor owns a piece of more
%          than one region
% 

%px:0   0   0   0   1   1   1   2   2   2   2   2   2   2   2   2   3   3   3   4   4   4   4   4   4
%px:                        2                                   3
%   |------ rx=0  ----------<------------- rx=1 ----------------<----------- rx=2 ------------------|            py  py
%   |                       |                                   |                                   |  ---
%   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24   |    0    0
%  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49   |    1    0
%  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74   |    2    0
%  75  76  77  78  79  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99 ry=0   3    0
% 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124   |    4    1
% 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149   |    5    1
% 150-151-152-153-154-155-156-157-158-159-160-161-162-163-164-165-166-167-168-169-170-171-172-173-174---<    6    1   2
% 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199   |    7    2
% 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224   |    8    2
% 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249   |    9    2
% 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274   |   10    2
% 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 ry=1  11    2
% 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324   |   12    2
% 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349   |   13    2
% 350 351 352 353 354 355 356 357 358 359 360 361 362 363 364 365 366 367 368 369 370 371 372 373 374   |   14    2
% 375-376-377-378-379-380-381-382-383-384-385-386-387-388-389-390-391-392-393-394-395-396-397-398-399---<   15    2   3
% 400 401 402 403 404 405 406 407 408 409 410 411 412 413 414 415 416 417 418 419 420 421 422 423 424   |   16    3
% 425 426 427 428 429 430 431 432 433 434 435 436 437 438 439 440 441 442 443 444 445 446 447 448 449   |   17    3
% 450 451 452 453 454 455 456 457 458 459 460 461 462 463 464 465 466 467 468 469 470 471 472 473 474   |   18    3
% 475 476 477 478 479 480 481 482 483 484 485 486 487 488 489 490 491 492 493 494 495 496 497 498 499 ry=2  19    4
% 500 501 502 503 504 505 506 507 508 509 510 511 512 513 514 515 516 517 518 519 520 521 522 523 524   |   20    4
% 525 526 527 528 529 530 531 532 533 534 535 536 537 538 539 540 541 542 543 544 545 546 547 548 549   |   21    4
% 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569 570 571 572 573 574   |   22    4
% 575 576 577 578 579 580 581 582 583 584 585 586 587 588 589 590 591 592 593 594 595 596 597 598 599   |   23    4
% 600 601 602 603 604 605 606 607 608 609 610 611 612 613 614 615 616 617 618 619 620 621 622 623 624   |   24    4
%   |                       |                                   |                                   |  ---
%   |-----------------------<-----------------------------------<-----------------------------------|
%
%   region = ry*3+rx  , proc = py*5+px
%
   whichCase = 'RegionsSpanProcs';
   rxInterfaceLocations = [6 15 24];   
   ryInterfaceLocations = [6 15 24];   
   pxProcChange         = [4  7 16 19 ];   % Note: the 7 and 16 (instead of 6 and 15)
   pyProcChange         = [4  7 16 19 ];   % at region interfaces where two procs share
                                           % a node. In general, we give the location
                                           % of the new processor's first non-shared point
elseif (strcmp(filename, 'caseSix') == true)  % caseSix

%px:0   0   0   0   0   0   0   0   0   1   1   1   1   1   2   2   2   2   3   3   3   3   3   3   3
%px:                                1                       1               2
%   |- rx=0 ----<-rx=1 -<- rx=2 ----<---- rx=3 -------------<-rx=4--<-rx=5--< ---- rx=6 ------------|            py  py
%   |           |       |           |                       |       |       |                       |     
%   0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24---<    0    0
%  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47  48  49 ry=0   1    0
%  50  51  52  53  54  55  56  57  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74   |    2    0
%--75--76--77--78--79--80--81--82--83--84--85--86--87--88--89--90--91--92--93--94--95--96--97--98--99---<    3    0
% 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 ry=1   4    0
%-125-126-127-128-129-130-131-132-133-134-135-136-137-138-139-140-141-142-143-144-145-146-147-148-149---<    5    0
% 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174   |    6    0    
% 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 ry=2   7    0
%-200-201-202-203-204-205-206-207-208-209-210-211-212-213-214-215-216-217-218-219-220-221-222-223-224---<    8    0   1
% 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249   |    9    1
% 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274   |   10    1
% 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 ry=3  11    1
% 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324   |   12    1
% 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349   |   13    1
%-350-351-352-353-354-355-356-357-358-359-360-361-362-363-364-365-366-367-368-369-370-371-372-373-374---<---14    2   1
% 375 376 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391 392 393 394 395 396 397 398 399 ry=4  15    2    
%-400-401-402-403-404-405-406-407-408-409-410-411-412-413-414-415-416-417-418-419-420-421-422-423-424---<---16    2
% 425 426 427 428 429 430 431 432 433 434 435 436 437 438 439 440 441 442 443 444 445 446 447 448 449 ry=5  17    2
%-450-451-452-453-454-455-456-457-458-459-460-461-462-463-464-465-466-467-468-469-470-471-472-473-474---<---18    3   2
% 475 476 477 478 479 480 481 482 483 484 485 486 487 488 489 490 491 492 493 494 495 496 497 498 499   |   19    3
% 500 501 502 503 504 505 506 507 508 509 510 511 512 513 514 515 516 517 518 519 520 521 522 523 524   |   20    3
% 525 526 527 528 529 530 531 532 533 534 535 536 537 538 539 540 541 542 543 544 545 546 547 548 549   |   21    3
% 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569 570 571 572 573 574 ry=6  22    3
% 575 576 577 578 579 580 581 582 583 584 585 586 587 588 589 590 591 592 593 594 595 596 597 598 599   |   23    3
% 600 601 602 603 604 605 606 607 608 609 610 611 612 613 614 615 616 617 618 619 620 621 622 623 624   |   24    3
%   |           |       |           |                       |       |       |                       |     
%   |- rx=0 ----<-rx=1 -<- rx=2 ----<---- rx=3 -------------<-rx=4--<-rx=5--< ---- rx=6 ------------|            
%px:0   0   0   0   0   0   0   0   0   1   1   1   1   1   2   2   2   2   3   3   3   3   3   3   3
%px:                                1                       1               2
%
%   region = ry*7+rx  , proc = py*4+px
%
   whichCase = 'MultipleRegionsPerProc';
   rxInterfaceLocations = [3 5 8 14 16 18 24];   
   ryInterfaceLocations = [3 5 8 14 16 18 24];   
   pxProcChange         = [9  15 19 ];   % Note: the 9 and 15 (instead of 8 and 14)
   pyProcChange         = [9  15 19 ];   % at region interfaces where two procs share
                                         % a node. In general, we give the location
                                         % of the new processor's first non-shared point

elseif(strcmp(filename, 'caseSeven') == true)  % caseSeven
  
%px:0   0   0   0   1   1   1   1
%px:            1                
%   |- rx=0 ----<-rx=1 -<- rx=2 |         py  py
%   |           |       |       |
%   0   1   2   3   4   5   6   7     0   0
%   8   9  10  11  12  13  14  15     1   0
%  16  17  18  19  20  21  22  23     2   0
%--24--25--26--27--28--29--30--31--<  3   0   1
%  32  33  34  35  36  37  38  39     4   1
%--40--41--42--43--44--45--46--47--<  5   1
%  48  49  50  51  52  53  54  55     6   1
%--56--57--58--59--60--61--62--63--<  7   1
  
   whichCase = 'MultipleRegionsPerProc';
   rxInterfaceLocations = [3 5 7];   
   ryInterfaceLocations = [3 5 7];   
   pxProcChange         = [3];
   pyProcChange         = [3];

elseif(strcmp(filename, 'caseEight') == true)  % caseEight
  
  % This is supposed to be the RegionsSpanProcs counterpart to caseSeven.
  error('Not implemented, yet.');
   
elseif(strcmp(filename, 'caseNine') == true)  % caseNine

%px:0   0   0   1   1
%px:        1                
%   |- rx=0 <-rx=1  |         py  py
%   |       |       |
%   0   1   2   3   4     0   0
%   5   6   7   8   9     1   0
%--10--11--12--13--14--<  2   0   1
%  15  16  17  18  19     3   1
%  20  21  22  23  24     4   1

  whichCase = 'MultipleRegionsPerProc';
  rxInterfaceLocations = [2 4];   
  ryInterfaceLocations = [2 4];   
  pxProcChange         = [3];
  pyProcChange         = [3];

elseif(strcmp(filename, 'caseFifteen') == true)  % caseFifteen

%px:0   0   0   0   1   1   1
%px:            1                
%   |- rx=0     <-rx=1      |   py  py
%   |           |           |
%   0   1   2   3   4   5   6   0   0
%   7   8   9  10  11  12  13  
%  14  15  16  17  18  19  20   2   0   1
%--21--22--23--24--25--26--27<   1
%  28  29  30  31  32  33  34   1
%  35  36  37  38  39  40  41
%  42  43  44  45  46  47  48

  whichCase = 'MultipleRegionsPerProc';
  rxInterfaceLocations = [3 6];   
  ryInterfaceLocations = [3 6];   
  pxProcChange         = [4];
  pyProcChange         = [4];
  
elseif(strcmp(filename, 'caseSixteen') == true)  % caseSixteen

%px:0   0   0   0   1   1   1
%px:            1                
%   |- rx=0     <-rx=1      |   py  py
%   |           |           |
%   0   1   2   3   4   5   6   0   0
%   7   8   9  10  11  12  13  
%  14  15  16  17  18  19  20   2   0   1
%--21--22--23--24--25--26--27<   1
%  28  29  30  31  32  33  34   1
%  35  36  37  38  39  40  41
%  42  43  44  45  46  47  48

  whichCase = 'MultipleRegionsPerProc';
  rxInterfaceLocations = [6 12];   
  ryInterfaceLocations = [6 12];   
  pxProcChange         = [7];
  pyProcChange         = [7];

elseif(strcmp(filename, 'caseSeventeen') == true)  % caseSeventeen

%px:0   0   0   0   1   1   1
%px:            1                
%   |- rx=0     <-rx=1      |   py  py
%   |           |           |
%   0   1   2   3   4   5   6   0   0
%   7   8   9  10  11  12  13  
%  14  15  16  17  18  19  20   2   0   1
%--21--22--23--24--25--26--27<   1
%  28  29  30  31  32  33  34   1
%  35  36  37  38  39  40  41
%  42  43  44  45  46  47  48

  whichCase = 'MultipleRegionsPerProc';
  rxInterfaceLocations = [30 60];   
  ryInterfaceLocations = [30 60];   
  pxProcChange         = [31];
  pyProcChange         = [31];
  
elseif(strcmp(filename, 'caseEightteen') == true)  % caseEightteen

%px:0   0   0   0   1   1   1
%px:            1                
%   |- rx=0     <-rx=1      |   py  py
%   |           |           |
%   0   1   2   3   4   5   6   0   0
%   7   8   9  10  11  12  13  
%  14  15  16  17  18  19  20   2   0   1
%--21--22--23--24--25--26--27<   1
%  28  29  30  31  32  33  34   1
%  35  36  37  38  39  40  41
%  42  43  44  45  46  47  48

  whichCase = 'MultipleRegionsPerProc';
  rxInterfaceLocations = [15 30 45];   
  ryInterfaceLocations = [15 30 45];   
  pxProcChange         = [16 31];
  pyProcChange         = [16 31];

elseif(strcmp(filename, 'caseNineteen') == true)  % caseNineteen

%px:0   0   0   0   1   1   1
%px:            1                
%   |- rx=0     <-rx=1      |   py  py
%   |           |           |
%   0   1   2   3   4   5   6   0   0
%   7   8   9  10  11  12  13  
%  14  15  16  17  18  19  20   2   0   1
%--21--22--23--24--25--26--27<   1
%  28  29  30  31  32  33  34   1
%  35  36  37  38  39  40  41
%  42  43  44  45  46  47  48

  whichCase = 'MultipleRegionsPerProc';
  rxInterfaceLocations = [15 30 45 60];   
  ryInterfaceLocations = [15 30 45 60];   
  pxProcChange         = [16 31 46];
  pyProcChange         = [16 31 46];
  
elseif(strcmp(filename, 'caseTwenty') == true)  % caseTwenty

%px:0   0   0   0   1   1   1
%px:            1                
%   |- rx=0     <-rx=1      |   py  py
%   |           |           |
%   0   1   2   3   4   5   6   0   0
%   7   8   9  10  11  12  13  
%  14  15  16  17  18  19  20   2   0   1
%--21--22--23--24--25--26--27<   1
%  28  29  30  31  32  33  34   1
%  35  36  37  38  39  40  41
%  42  43  44  45  46  47  48

  whichCase = 'MultipleRegionsPerProc';
  rxInterfaceLocations = [9 18];   
  ryInterfaceLocations = [9 18];   
  pxProcChange         = [10];
  pyProcChange         = [10];

end

nrx = length(rxInterfaceLocations);
nry = length(ryInterfaceLocations);
nRegions = nrx*nry;
nx       = rxInterfaceLocations(nrx)+1;
ny       = ryInterfaceLocations(nry)+1;
nNodes   = nx*ny;
npx      = length(pxProcChange)+1;
npy      = length(pyProcChange)+1;
nProcs   = npx*npy;
globalDims=  -ones(nRegions,2);
procsCorner = -ones(nRegions,nProcs,2);
relCorner   = -ones(nRegions,nProcs,2);
regionCorner= -ones(nRegions,2);
localDims =  -ones(nRegions,nProcs,2);

fp = fopen(filename,'w');
if fp ~= -1, 
    fprintf(fp,'     #nodes  #regions   #procs   #which case\n');
    fprintf(fp,'%8d %8d %8d       %s\n',nNodes,nRegions,nProcs,whichCase);
    fprintf(fp,'     nodeID   regionID   procID\n');

    % basic idea here is to loop over every grid point while at the same time we keep
    % track of which proc is assigned to the grid point via the code lines 
    %
    %   if (pj ~= npy-1)&&(j == pyProcChange(pj+1)), pj = pj+1; end;
    %   if (pi ~= npx-1)&&(i == pxProcChange(pi+1)), pi = pi+1; end;
    %
    % as well as tracking which region is associated with via code lines
    %
    %   if j == ryInterfaceLocations(rj+1)+1, rj = rj+1; end;
    %   if i == rxInterfaceLocations(ri+1)+1, ri = ri+1; end;
    %
    % Then, we print/compute box info for the grid point. When a grid point is
    % on the edge of a region interface, then it lives in multiple regions and 
    % might live on multiple procs. So the if's that follows first fprintf/computeBoxDims()
    % are intended to handle these region interface cases.

    pj = 0;  rj = 0;
    for j=0:ny-1, 
       pi = 0; ri = 0;
       if (pj ~= npy-1)&&(j == pyProcChange(pj+1)), pj = pj+1; end;
       if j == ryInterfaceLocations(rj+1)+1, rj = rj+1; end;
               % at interfaces, we defer the region switch to the 
               % first non-shared node, hence the +1 above.
       for i=0:nx-1,
          if (pi ~= npx-1)&&(i == pxProcChange(pi+1)), pi = pi+1; end;
          if i == rxInterfaceLocations(ri+1)+1, ri = ri+1; end;
          me     =  j*nx +  i;
          proc   = pj*npx + pi;
          region = rj*nrx + ri;
          fprintf(fp,'%8d  %8d %8d\n',me,region,proc); %,pi,pj,ri,rj);
          [procsCorner,globalDims,localDims,regionCorner] = computeBoxDims(region+1,proc+1,i,j,procsCorner,globalDims,localDims,regionCorner);
          tempri = ri; temprj = rj; temppi = pi; temppj = pj;
          if (i == rxInterfaceLocations(ri+1)) && (ri ~= nrx-1),
             tempri = ri+1;
             if (pi ~= npx-1)&&(i == pxProcChange(pi+1)-1), temppi = pi+1; end;
             proc   = temppj*npx + temppi;
             region = temprj*nrx + tempri;
             fprintf(fp,'%8d  %8d %8d\n',me,region,proc); %,temppi,temppj,tempri,temprj);
             [procsCorner,globalDims,localDims,regionCorner] = computeBoxDims(region+1,proc+1,i,j,procsCorner,globalDims,localDims,regionCorner);
          end;
          tempri = ri; temprj = rj; temppi = pi; temppj = pj;
          if j == ryInterfaceLocations(rj+1) && (rj ~= nry-1),
             temprj = rj+1;
             if (pj ~= npy-1)&&(j == pyProcChange(pj+1)-1), temppj = pj+1; end;
             proc   = temppj*npx + temppi;
             region = temprj*nrx + tempri;
             fprintf(fp,'%8d  %8d %8d\n',me,region,proc); %,temppi,temppj,tempri,temprj);
             [procsCorner,globalDims,localDims,regionCorner] = computeBoxDims(region+1,proc+1,i,j,procsCorner,globalDims,localDims,regionCorner);
          end;
          tempri = ri; temprj = rj; temppi = pi; temppj = pj;
          if (i == rxInterfaceLocations(ri+1)) && (j == ryInterfaceLocations(rj+1)) ...
              && (ri ~= nrx-1) && (rj ~= nry-1),
             tempri = ri+1; temprj = rj+1;
             if (pi ~= npx-1)&&(i == pxProcChange(pi+1)-1), temppi = pi+1; end;
             if (pj ~= npy-1)&&(j == pyProcChange(pj+1)-1), temppj = pj+1; end;
             proc   = temppj*npx + temppi;
             region = temprj*nrx + tempri;
             fprintf(fp,'%8d  %8d %8d\n',me,region,proc); %,temppi,temppj,tempri,temprj);
             [procsCorner,globalDims,localDims,regionCorner] = computeBoxDims(region+1,proc+1,i,j,procsCorner,globalDims,localDims,regionCorner);
          end
       end;
    end; 
    fclose(fp);
else
    fprintf('could not open file %s\n',filename);
    keyboard;
end;

%
%  now subtract the corner position associated with a proc/region from the true corner position
%  of a region (over all procs).
%
for i=1:nProcs,
    inds = find(procsCorner(:,i,1) ~= -1);
    relCorner(inds,i,1) = procsCorner(inds,i,1) - regionCorner(inds,1);
    relCorner(inds,i,2) = procsCorner(inds,i,2) - regionCorner(inds,2);
end

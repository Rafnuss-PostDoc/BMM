Track table
| Field Name	|  Description	| 
| --- 			|---|
| Code 			| according to table Site |
| Date 			| Date and time |
| FNr 			| Flight path number (starts for a day at 1) |
| z 			| Average altitude (m) |
| Rg 			| Mean flight direction (°) |
| Ra 			| Mean air direction (°) |
| Vg 			| Mean groundspeed (cm / s) |
| Va 			| Mean airspeed (cm / s) |
| Vz 			| Mean rate of climb (cm / s) |
| lnt 			| Number of recorded intervals |
| MVr 			| Average effective speed (straightness, 1000 = absolutely straight) |
| mRab 			| Average tendency to turn per 20s (') |
| number 		| Number code (O = undetermined / wind; 1 = 1; 2 = 1-2; 3 = 3-5; 4 = 6-20; 5 = 21-50; 6 = 51-100; 7> 100} |
| FieldClass 	| Wing flapping class {l = large water bird {<9Hz); 2 = small water bird; 3 = large songbird {<12Hz); 4 = small songbird; 5 = common swift; 6 = bird of prey; 7 = big bird {<6Hz); 9 = indefinite / wind) |
| Species 		| Species code |
| F_field 		| Field wing beat frequency (Hz / 10) |
| Rw 			| Mean wind direction (°) |
| Vw 			| Mean wind speed (cm / s) |
| MDVw 			| Mean deviation of the rate of climb of the wind balloon (cm / s) |
| Press 		| Pressure |
| Temp 			| Temperature |
| Humid 		| Moisture |
| DW 			| Mean horizontal distance for wind measurement (m) |
| MDzW 			| Mean vertical distance for wind measurement (m) |
| MtDW 			| Time difference to wind measurement (hhmm) |
| MtDS 			| Temporal difference to the probe (hhmm) |
| AltitudeSite 	| Height of the observation point in meters above sea level |

Site table
| Field Name 		| Description |
| --- 				| --- |
| Code 				| Project code |
| SiteID 			| Unique number for a series of measurements at one location; Used for record numbering in the radar databases |
| RadarID 			| Radar device used |
| Name 				| Location name |
| Longitude 		| Longitude in 'East |
| Latitude 			| Latitude in 'North |
| Altitude 			| Height of the observation site in m above sea level. M. |
| ProjectStart		| Date when field data registration started |
| ProjectEnd 		| Date at the end of field data registration |
| Timezone 			| Deviation of the time from UTC in h |
| ResponsiblePerson | Person responsible for field work |
| Sitedescription 	| Description of the location |
| Remarks 			| |

Wind table
| Field Name 	| Description |
| --- 			| --- |
| Code 			| Location code |
| DateT ime 	| Date and time at the start of the measurement |
| lnt 			| Sequence number of the 20s interval for averaging the registered second data |
| z 			| Height above location in meters |
| Dg 			| Direction of flight of the anemometer balloon in '(> wind direction -180') |
| Vg 			| Airspeed of the anemometer balloon in cm / s (> wind speed) |
| Vz 			| Vertical speed of the anemometer balloon in cm / s (> rate of climb) |
| Points 		| Number of second data taken into account for the calculation of the mean values ​​in this interval (maximum 21; If the speed calculation for the last second resulted in more than SOm / s, the second values ​​were not recognized as valid; ROHDAT!) |
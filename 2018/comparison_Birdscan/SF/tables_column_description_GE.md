Track table
| Field Name	| Description	| 
| --- 			|---|
| Code			| gemäss Tabelle Site |
| Date			| Datum und Uhrzeit |
| FNr			| Flugwegnummer (beginnt jeweils tageweise bei 1} |
| z				| Mittlere Flughöhe (m) |
| Rg			| Mittlere Flugrichtung (') |
| Ra			| Mittlere Eigenrichtung (') |
| Vg			| Mittlere Fluggeschwindigkeit (cm/s) |
| Va			| Mittlere Eigengeschwindigkeit (cm/s) |
| Vz			| Mittlere Steigrate (cm/s) |
| lnt			| Anzahl aufgezeicnete Intervalle |
| MVr			| Mittlere effektive Geschwindigkeit (Geradlinigkeit, 1000 = absolut geradlinig) |
| mRab			| Mittlere Abbiegetendenz pro 20s (') |
| number		| Anzahlcode (O=unbestimmt/Wind; 1=1; 2=1-2; 3=3-5; 4=6-20; 5=21-50; 6=51-100; 7>100} |
| FieldClass	| Flügelschlagklasse {l=grosser Wasservogel {<9Hz); 2=kleiner Wasservogel; 3=grosser Singvogel {<12Hz); 4=kleiner Singvogel; 5=Mauersegler; 6=Greifvogel; 7=grosser Vogel {<6Hz); 9=unbestimmt/Wind) |
| Species		| Artcode |
| F_field		| Feld-Flügelschlagfrequenz (Hz/10) |
| Rw			| Mittlere Windrichtung (') |
| Vw			| Mittlere Windgeschwindigkeit (cm/s) |
| MDVw			| Mittlere Abweichung der Steiggeschwindigkeit des Windballons (cm/s) |
| Press			| Druck |
| Temp			| Temperatur |
| Humid			| Feuchtigkeit |
| DW			| Mittlere Horizontaldist anz zur Windmessung (m) |
| MDzW			| Mittlere Vertikaldist anz zur Winmessung (m) |
| MtDW			| Ze itliche Differenz zur Windmessung (hhmm) |
| MtDS			| Zeitliche Differenz zur Sondage (hhmm) |
| AltitudeSite	| Höhe des Beobachtungsortes in Meter über Meer |

Site table
| Field Name		| Description	| 
| --- 				|---|
| Code				| Projectcode |
| SitelD			| Eindeutige Nummer für eine Messreihe an einem Ort; Wird gebraucht für Recordnummerierung in den Radardatenbanken |
| RadarlD			| Verwendetes Radargerät |
| Name				| Standortbezeichnung |
| Longitude			| Geografische Länge in ' Ost |
| Latitude			| Geografische Breite in 'Nord |
| Altitude			| Höhe des Beobachtungsortes in m ü. M. |
| Proj ectStart		| Datum bei Start der Felddatenregistrierung |
| ProjectEnd		| Datum bei Ende der Felddatenregistr ierung |
| Timezone			| Abweichung der Zeitangabe gegenüber UTC in h |
| ResponsiblePerson	| Für die Feldarbeit verantwortliche Person |
| Sitedescription	| Beschre ibung des Standortes |
| Remarks			|  |

Wind table
| Field Name	| Description	| 
| --- 			|---|
| Code			| Standortcode |
| DateT ime		| Datum und Zeit bei Beginn der Messung |
| lnt			| Laufnummer des 20s-lntervalls für Mittelung der registrierten Sekundendaten |
| z				| Höhe über Standort in Meter |
| Dg			| Flugrichtung des Windmessballons in '(> Windrichtung -180') |
| Vg			| Fluggeschwindigkeit des Windmessballons in cm/s (> Windgeschwindigkeit) |
| Vz			| Vertikalgeschwindigkeit des Windmessballons in cm/s (>Steigrate) |
| Points		| Anzahl für die Berechnung der Mittelwerte berücksichtigte Sekundendaten in diesem Intervall (maximal 21; Wenn die Geschwindigkeitesberechnung für die jeweils letzte Sekunde mehr alsSOm/ s ergab wurden die Sekundenwerte nicht als gültig anerkannt; ROHDAT!) |

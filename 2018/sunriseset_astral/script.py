import datetime
import astral
import json

# latlon = [[51.1917, 3.0642],[49.9143, 5.5056],[49.6583, 13.8178],[49.5011, 16.7885],[53.564, 6.7483],[54.0044, 10.0468],[51.1246, 13.7686],[49.5407, 12.4028],[53.3394, 7.025],[51.4055, 6.9669],[51.3112, 8.802],[47.8736, 8.0036],[52.4601, 9.6945],[48.0431, 10.2204],[50.5001, 11.135],[50.1097, 6.5485],[49.9859, 8.714],[52.6486, 13.8578],[54.1757, 12.0581],[48.1747, 12.1018],[48.5853, 9.7828],[52.1601, 11.1761],[55.11275, 14.887517],[55.1731110, 8.552],[57.4893060, 10.136472],[55.3261940, 12.449278],[56.0240069, 10.02459060],[59.3976700119674, 24.6021000109613],[58.4823100268841, 25.5186601169407],[36.8325000442565, -2.08221998065711],[39.4288901053369, -6.28527989611031],[41.408060118556, 1.8847200460732],[43.1688901036978, -8.52693993598224],[41.9955601170659, -4.60277996957303],[28.0186100117862, -15.6144398823381],[40.1758299954236, -3.71360991150142],[36.613330040127, -4.65916989371182],[38.264440074563, -1.18971996009351],[39.379720017314, 2.78499998152256],[43.4625001251698, -6.30193993449212],[37.6875001005828, -6.33443992584945],[43.403330091387, -2.84193992614748],[39.1761100105941, -0.252109877765192],[41.7338900081813, -0.545830037444844],[60.9038700163364, 27.1080600656569],[61.7673399858177, 23.0764499679208],[61.9069920480251, 29.7977200895548],[60.1284700445831, 21.6433800570667],[62.8626000881195, 27.3814700357616],[67.1391000598669, 26.8969160690904],[64.7749301232398, 26.3188800774515],[60.2706200815737, 24.8690200410783],[63.1048401072621, 23.8208600319922],[50.13583, 1.83472],[42.12972, 9.49639],[50.12833, 3.81194],[47.35528, 4.77583],[44.32306, 4.76222],[44.83139, -0.69167],[47.05861, 2.35944],[48.92722, -0.14944],[46.69861, 0.06556],[43.21667, 6.37278],[45.10444, 1.36944],[45.29, 3.70944],[43.99056, 2.60972],[43.62472, -0.60917],[47.36861, 7.01917],[48.71583, 6.58167],[43.80611, 4.50278],[46.06778, 4.44528],[42.91833, 2.865],[48.46083, -4.43],[43.57444, 1.37611],[48.77389, 2.0075],[47.3375, -1.65639],[48.46222, 4.30944],[45.8834500797093, 17.2008900716901],[45.5027001164854, 18.5613001137972],[52.95279, 4.79061],[51.8369, 5.1381],[50.39417, 20.07972],[54.38425, 18.45631],[52.4052190, 20.960911],[50.892, 16.0395],[52.41326, 16.79706],[50.15167, 18.72667],[50.11409, 22.03704],[53.79028, 15.83111],[56.3675003051758, 12.8516998291016],[59.6543998718, 17.9463005066],[57.3033981323242, 18.4002990722656],[61.5771102905273, 16.7144336700439],[67.7088012695, 20.6177997589],[56.295501709, 15.6103000641],[60.7230415344238, 14.877571105957],[65.4309005737, 21.8649997711],[63.6394691467285, 18.4018821716309],[63.2949981689453, 14.7590999603272],[58.2555999755859, 12.826000213623],[58.1058998108, 15.9363002777],[46.0677699744701, 15.284899994731],[46.097980029881, 14.2282400652766],[48.2561, 17.1531],[48.7828950, 20.987277],[49.27167, 19.24935],[48.24043, 19.25743]]

# t0 = datetime.date(2018, 1, 1)

# d=[]
# for ll in latlon:
#     d2=[]
#     a = astral.Location(('name', 'region', ll[0], ll[1], 'UTC', 0))
#     for t in range(0,366):
#         try:
#             tmp = a.sun(t0 + datetime.timedelta(t))
#             d2.append([str(tmp['dawn']).split('+')[0], str(tmp['sunrise']).split('+')[0], str(tmp['sunset']).split('+')[0], str(tmp['dusk']).split('+')[0]])
#         except:
#             d2.append([])
#     d.append(d2)

# a.twilight(t0+datetime.timedelta(t))
# a.dusk(t0 + datetime.timedelta(t))
# a.sunset(t0 + datetime.timedelta(t))
# a.sunrise(t0 + datetime.timedelta(t))


# with open('sunriseset_latlon.json', 'w') as fout:
#     json.dump(d, fout)





latlon = [[48.5,-4.75],[48,-4.5],[48.25,-4.5],[48.5,-4.5],[48,-4.25],[48.25,-4.25],[48.5,-4.25],[48,-4],[48.25,-4],[48.5,-4],[48,-3.75],[48.25,-3.75],[48.5,-3.75],[47.75,-3.5],[48,-3.5],[48.25,-3.5],[48.5,-3.5],[48.75,-3.5],[47.75,-3.25],[48,-3.25],[48.25,-3.25],[48.5,-3.25],[48.75,-3.25],[47.75,-3],[48,-3],[48.25,-3],[48.5,-3],[48.75,-3],[47.75,-2.75],[48,-2.75],[48.25,-2.75],[48.5,-2.75],[47.5,-2.5],[47.75,-2.5],[48,-2.5],[48.25,-2.5],[48.5,-2.5],[43,-2.25],[43.25,-2.25],[47.5,-2.25],[47.75,-2.25],[48,-2.25],[48.25,-2.25],[48.5,-2.25],[43,-2],[43.25,-2],[46.75,-2],[47,-2],[47.25,-2],[47.5,-2],[47.75,-2],[48,-2],[48.25,-2],[48.5,-2],[43,-1.75],[43.25,-1.75],[46.5,-1.75],[46.75,-1.75],[47,-1.75],[47.25,-1.75],[47.5,-1.75],[47.75,-1.75],[48,-1.75],[48.25,-1.75],[48.5,-1.75],[49.5,-1.75],[43,-1.5],[43.25,-1.5],[43.5,-1.5],[46.5,-1.5],[46.75,-1.5],[47,-1.5],[47.25,-1.5],[47.5,-1.5],[47.75,-1.5],[48,-1.5],[48.25,-1.5],[48.5,-1.5],[48.75,-1.5],[49,-1.5],[49.25,-1.5],[49.5,-1.5],[43,-1.25],[43.25,-1.25],[43.5,-1.25],[43.75,-1.25],[44,-1.25],[44.25,-1.25],[46.25,-1.25],[46.5,-1.25],[46.75,-1.25],[47,-1.25],[47.25,-1.25],[47.5,-1.25],[47.75,-1.25],[48,-1.25],[48.25,-1.25],[48.5,-1.25],[48.75,-1.25],[49,-1.25],[49.25,-1.25],[49.5,-1.25],[43,-1],[43.25,-1],[43.5,-1],[43.75,-1],[44,-1],[44.25,-1],[44.5,-1],[44.75,-1],[45,-1],[45.25,-1],[45.5,-1],[45.75,-1],[46,-1],[46.25,-1],[46.5,-1],[46.75,-1],[47,-1],[47.25,-1],[47.5,-1],[47.75,-1],[48,-1],[48.25,-1],[48.5,-1],[48.75,-1],[49,-1],[49.25,-1],[43,-0.75],[43.25,-0.75],[43.5,-0.75],[43.75,-0.75],[44,-0.75],[44.25,-0.75],[44.5,-0.75],[44.75,-0.75],[45,-0.75],[45.25,-0.75],[45.5,-0.75],[45.75,-0.75],[46,-0.75],[46.25,-0.75],[46.5,-0.75],[46.75,-0.75],[47,-0.75],[47.25,-0.75],[47.5,-0.75],[47.75,-0.75],[48,-0.75],[48.25,-0.75],[48.5,-0.75],[48.75,-0.75],[49,-0.75],[49.25,-0.75],[43,-0.5],[43.25,-0.5],[43.5,-0.5],[43.75,-0.5],[44,-0.5],[44.25,-0.5],[44.5,-0.5],[44.75,-0.5],[45,-0.5],[45.25,-0.5],[45.5,-0.5],[45.75,-0.5],[46,-0.5],[46.25,-0.5],[46.5,-0.5],[46.75,-0.5],[47,-0.5],[47.25,-0.5],[47.5,-0.5],[47.75,-0.5],[48,-0.5],[48.25,-0.5],[48.5,-0.5],[48.75,-0.5],[49,-0.5],[49.25,-0.5],[43,-0.25],[43.25,-0.25],[43.5,-0.25],[43.75,-0.25],[44,-0.25],[44.25,-0.25],[44.5,-0.25],[44.75,-0.25],[45,-0.25],[45.25,-0.25],[45.5,-0.25],[45.75,-0.25],[46,-0.25],[46.25,-0.25],[46.5,-0.25],[46.75,-0.25],[47,-0.25],[47.25,-0.25],[47.5,-0.25],[47.75,-0.25],[48,-0.25],[48.25,-0.25],[48.5,-0.25],[48.75,-0.25],[49,-0.25],[49.25,-0.25],[43,0],[43.25,0],[43.5,0],[43.75,0],[44,0],[44.25,0],[44.5,0],[44.75,0],[45,0],[45.25,0],[45.5,0],[45.75,0],[46,0],[46.25,0],[46.5,0],[46.75,0],[47,0],[47.25,0],[47.5,0],[47.75,0],[48,0],[48.25,0],[48.5,0],[48.75,0],[49,0],[49.25,0],[43,0.25],[43.25,0.25],[43.5,0.25],[43.75,0.25],[44,0.25],[44.25,0.25],[44.5,0.25],[44.75,0.25],[45,0.25],[45.25,0.25],[45.5,0.25],[45.75,0.25],[46,0.25],[46.25,0.25],[46.5,0.25],[46.75,0.25],[47,0.25],[47.25,0.25],[47.5,0.25],[47.75,0.25],[48,0.25],[48.25,0.25],[48.5,0.25],[48.75,0.25],[49,0.25],[49.25,0.25],[49.5,0.25],[43,0.5],[43.25,0.5],[43.5,0.5],[43.75,0.5],[44,0.5],[44.25,0.5],[44.5,0.5],[44.75,0.5],[45,0.5],[45.25,0.5],[45.5,0.5],[45.75,0.5],[46,0.5],[46.25,0.5],[46.5,0.5],[46.75,0.5],[47,0.5],[47.25,0.5],[47.5,0.5],[47.75,0.5],[48,0.5],[48.25,0.5],[48.5,0.5],[48.75,0.5],[49,0.5],[49.25,0.5],[49.5,0.5],[49.75,0.5],[43,0.75],[43.25,0.75],[43.5,0.75],[43.75,0.75],[44,0.75],[44.25,0.75],[44.5,0.75],[44.75,0.75],[45,0.75],[45.25,0.75],[45.5,0.75],[45.75,0.75],[46,0.75],[46.25,0.75],[46.5,0.75],[46.75,0.75],[47,0.75],[47.25,0.75],[47.5,0.75],[47.75,0.75],[48,0.75],[48.25,0.75],[48.5,0.75],[48.75,0.75],[49,0.75],[49.25,0.75],[49.5,0.75],[49.75,0.75],[43,1],[43.25,1],[43.5,1],[43.75,1],[44,1],[44.25,1],[44.5,1],[44.75,1],[45,1],[45.25,1],[45.5,1],[45.75,1],[46,1],[46.25,1],[46.5,1],[46.75,1],[47,1],[47.25,1],[47.5,1],[47.75,1],[48,1],[48.25,1],[48.5,1],[48.75,1],[49,1],[49.25,1],[49.5,1],[49.75,1],[43,1.25],[43.25,1.25],[43.5,1.25],[43.75,1.25],[44,1.25],[44.25,1.25],[44.5,1.25],[44.75,1.25],[45,1.25],[45.25,1.25],[45.5,1.25],[45.75,1.25],[46,1.25],[46.25,1.25],[46.5,1.25],[46.75,1.25],[47,1.25],[47.25,1.25],[47.5,1.25],[47.75,1.25],[48,1.25],[48.25,1.25],[48.5,1.25],[48.75,1.25],[49,1.25],[49.25,1.25],[49.5,1.25],[49.75,1.25],[43,1.5],[43.25,1.5],[43.5,1.5],[43.75,1.5],[44,1.5],[44.25,1.5],[44.5,1.5],[44.75,1.5],[45,1.5],[45.25,1.5],[45.5,1.5],[45.75,1.5],[46,1.5],[46.25,1.5],[46.5,1.5],[46.75,1.5],[47,1.5],[47.25,1.5],[47.5,1.5],[47.75,1.5],[48,1.5],[48.25,1.5],[48.5,1.5],[48.75,1.5],[49,1.5],[49.25,1.5],[49.5,1.5],[49.75,1.5],[50,1.5],[43,1.75],[43.25,1.75],[43.5,1.75],[43.75,1.75],[44,1.75],[44.25,1.75],[44.5,1.75],[44.75,1.75],[45,1.75],[45.25,1.75],[45.5,1.75],[45.75,1.75],[46,1.75],[46.25,1.75],[46.5,1.75],[46.75,1.75],[47,1.75],[47.25,1.75],[47.5,1.75],[47.75,1.75],[48,1.75],[48.25,1.75],[48.5,1.75],[48.75,1.75],[49,1.75],[49.25,1.75],[49.5,1.75],[49.75,1.75],[50,1.75],[50.25,1.75],[50.5,1.75],[50.75,1.75],[43,2],[43.25,2],[43.5,2],[43.75,2],[44,2],[44.25,2],[44.5,2],[44.75,2],[45,2],[45.25,2],[45.5,2],[45.75,2],[46,2],[46.25,2],[46.5,2],[46.75,2],[47,2],[47.25,2],[47.5,2],[47.75,2],[48,2],[48.25,2],[48.5,2],[48.75,2],[49,2],[49.25,2],[49.5,2],[49.75,2],[50,2],[50.25,2],[50.5,2],[50.75,2],[43,2.25],[43.25,2.25],[43.5,2.25],[43.75,2.25],[44,2.25],[44.25,2.25],[44.5,2.25],[44.75,2.25],[45,2.25],[45.25,2.25],[45.5,2.25],[45.75,2.25],[46,2.25],[46.25,2.25],[46.5,2.25],[46.75,2.25],[47,2.25],[47.25,2.25],[47.5,2.25],[47.75,2.25],[48,2.25],[48.25,2.25],[48.5,2.25],[48.75,2.25],[49,2.25],[49.25,2.25],[49.5,2.25],[49.75,2.25],[50,2.25],[50.25,2.25],[50.5,2.25],[50.75,2.25],[43,2.5],[43.25,2.5],[43.5,2.5],[43.75,2.5],[44,2.5],[44.25,2.5],[44.5,2.5],[44.75,2.5],[45,2.5],[45.25,2.5],[45.5,2.5],[45.75,2.5],[46,2.5],[46.25,2.5],[46.5,2.5],[46.75,2.5],[47,2.5],[47.25,2.5],[47.5,2.5],[47.75,2.5],[48,2.5],[48.25,2.5],[48.5,2.5],[48.75,2.5],[49,2.5],[49.25,2.5],[49.5,2.5],[49.75,2.5],[50,2.5],[50.25,2.5],[50.5,2.5],[50.75,2.5],[51,2.5],[43,2.75],[43.25,2.75],[43.5,2.75],[43.75,2.75],[44,2.75],[44.25,2.75],[44.5,2.75],[44.75,2.75],[45,2.75],[45.25,2.75],[45.5,2.75],[45.75,2.75],[46,2.75],[46.25,2.75],[46.5,2.75],[46.75,2.75],[47,2.75],[47.25,2.75],[47.5,2.75],[47.75,2.75],[48,2.75],[48.25,2.75],[48.5,2.75],[48.75,2.75],[49,2.75],[49.25,2.75],[49.5,2.75],[49.75,2.75],[50,2.75],[50.25,2.75],[50.5,2.75],[50.75,2.75],[51,2.75],[43,3],[43.25,3],[43.5,3],[43.75,3],[44,3],[44.25,3],[44.5,3],[44.75,3],[45,3],[45.25,3],[45.5,3],[45.75,3],[46,3],[46.25,3],[46.5,3],[46.75,3],[47,3],[47.25,3],[47.5,3],[47.75,3],[48,3],[48.25,3],[48.5,3],[48.75,3],[49,3],[49.25,3],[49.5,3],[49.75,3],[50,3],[50.25,3],[50.5,3],[50.75,3],[51,3],[51.25,3],[43.25,3.25],[43.5,3.25],[43.75,3.25],[44,3.25],[44.25,3.25],[44.5,3.25],[44.75,3.25],[45,3.25],[45.25,3.25],[45.5,3.25],[45.75,3.25],[46,3.25],[46.25,3.25],[46.5,3.25],[46.75,3.25],[47,3.25],[47.25,3.25],[47.5,3.25],[47.75,3.25],[48,3.25],[48.25,3.25],[48.5,3.25],[48.75,3.25],[49,3.25],[49.25,3.25],[49.5,3.25],[49.75,3.25],[50,3.25],[50.25,3.25],[50.5,3.25],[50.75,3.25],[51,3.25],[51.25,3.25],[43.25,3.5],[43.5,3.5],[43.75,3.5],[44,3.5],[44.25,3.5],[44.5,3.5],[44.75,3.5],[45,3.5],[45.25,3.5],[45.5,3.5],[45.75,3.5],[46,3.5],[46.25,3.5],[46.5,3.5],[46.75,3.5],[47,3.5],[47.25,3.5],[47.5,3.5],[47.75,3.5],[48,3.5],[48.25,3.5],[48.5,3.5],[48.75,3.5],[49,3.5],[49.25,3.5],[49.5,3.5],[49.75,3.5],[50,3.5],[50.25,3.5],[50.5,3.5],[50.75,3.5],[51,3.5],[51.25,3.5],[43.5,3.75],[43.75,3.75],[44,3.75],[44.25,3.75],[44.5,3.75],[44.75,3.75],[45,3.75],[45.25,3.75],[45.5,3.75],[45.75,3.75],[46,3.75],[46.25,3.75],[46.5,3.75],[46.75,3.75],[47,3.75],[47.25,3.75],[47.5,3.75],[47.75,3.75],[48,3.75],[48.25,3.75],[48.5,3.75],[48.75,3.75],[49,3.75],[49.25,3.75],[49.5,3.75],[49.75,3.75],[50,3.75],[50.25,3.75],[50.5,3.75],[50.75,3.75],[51,3.75],[51.25,3.75],[51.5,3.75],[43.5,4],[43.75,4],[44,4],[44.25,4],[44.5,4],[44.75,4],[45,4],[45.25,4],[45.5,4],[45.75,4],[46,4],[46.25,4],[46.5,4],[46.75,4],[47,4],[47.25,4],[47.5,4],[47.75,4],[48,4],[48.25,4],[48.5,4],[48.75,4],[49,4],[49.25,4],[49.5,4],[49.75,4],[50,4],[50.25,4],[50.5,4],[50.75,4],[51,4],[51.25,4],[51.5,4],[51.75,4],[43.5,4.25],[43.75,4.25],[44,4.25],[44.25,4.25],[44.5,4.25],[44.75,4.25],[45,4.25],[45.25,4.25],[45.5,4.25],[45.75,4.25],[46,4.25],[46.25,4.25],[46.5,4.25],[46.75,4.25],[47,4.25],[47.25,4.25],[47.5,4.25],[47.75,4.25],[48,4.25],[48.25,4.25],[48.5,4.25],[48.75,4.25],[49,4.25],[49.25,4.25],[49.5,4.25],[49.75,4.25],[50,4.25],[50.25,4.25],[50.5,4.25],[50.75,4.25],[51,4.25],[51.25,4.25],[51.5,4.25],[51.75,4.25],[52,4.25],[44,4.5],[44.25,4.5],[44.5,4.5],[44.75,4.5],[45,4.5],[45.25,4.5],[45.5,4.5],[45.75,4.5],[46,4.5],[46.25,4.5],[46.5,4.5],[46.75,4.5],[47,4.5],[47.25,4.5],[47.5,4.5],[47.75,4.5],[48,4.5],[48.25,4.5],[48.5,4.5],[48.75,4.5],[49,4.5],[49.25,4.5],[49.5,4.5],[49.75,4.5],[50,4.5],[50.25,4.5],[50.5,4.5],[50.75,4.5],[51,4.5],[51.25,4.5],[51.5,4.5],[51.75,4.5],[52,4.5],[52.25,4.5],[44.25,4.75],[44.5,4.75],[44.75,4.75],[45,4.75],[45.25,4.75],[45.5,4.75],[45.75,4.75],[46,4.75],[46.25,4.75],[46.5,4.75],[46.75,4.75],[47,4.75],[47.25,4.75],[47.5,4.75],[47.75,4.75],[48,4.75],[48.25,4.75],[48.5,4.75],[48.75,4.75],[49,4.75],[49.25,4.75],[49.5,4.75],[49.75,4.75],[50,4.75],[50.25,4.75],[50.5,4.75],[50.75,4.75],[51,4.75],[51.25,4.75],[51.5,4.75],[51.75,4.75],[52,4.75],[52.25,4.75],[52.5,4.75],[52.75,4.75],[44.5,5],[44.75,5],[45,5],[45.25,5],[45.5,5],[45.75,5],[46,5],[46.25,5],[46.5,5],[46.75,5],[47,5],[47.25,5],[47.5,5],[47.75,5],[48,5],[48.25,5],[48.5,5],[48.75,5],[49,5],[49.25,5],[49.5,5],[49.75,5],[50,5],[50.25,5],[50.5,5],[50.75,5],[51,5],[51.25,5],[51.5,5],[51.75,5],[52,5],[52.25,5],[52.5,5],[52.75,5],[44.75,5.25],[45,5.25],[45.25,5.25],[45.5,5.25],[45.75,5.25],[46,5.25],[46.25,5.25],[46.5,5.25],[46.75,5.25],[47,5.25],[47.25,5.25],[47.5,5.25],[47.75,5.25],[48,5.25],[48.25,5.25],[48.5,5.25],[48.75,5.25],[49,5.25],[49.25,5.25],[49.5,5.25],[49.75,5.25],[50,5.25],[50.25,5.25],[50.5,5.25],[50.75,5.25],[51,5.25],[51.25,5.25],[51.5,5.25],[51.75,5.25],[52,5.25],[52.25,5.25],[52.5,5.25],[52.75,5.25],[45,5.5],[45.25,5.5],[45.5,5.5],[45.75,5.5],[46,5.5],[46.25,5.5],[46.5,5.5],[46.75,5.5],[47,5.5],[47.25,5.5],[47.5,5.5],[47.75,5.5],[48,5.5],[48.25,5.5],[48.5,5.5],[48.75,5.5],[49,5.5],[49.25,5.5],[49.5,5.5],[49.75,5.5],[50,5.5],[50.25,5.5],[50.5,5.5],[50.75,5.5],[51,5.5],[51.25,5.5],[51.5,5.5],[51.75,5.5],[52,5.5],[52.25,5.5],[52.5,5.5],[52.75,5.5],[53,5.5],[53.25,5.5],[45.25,5.75],[45.5,5.75],[45.75,5.75],[46,5.75],[46.25,5.75],[46.5,5.75],[46.75,5.75],[47,5.75],[47.25,5.75],[47.5,5.75],[47.75,5.75],[48,5.75],[48.25,5.75],[48.5,5.75],[48.75,5.75],[49,5.75],[49.25,5.75],[49.5,5.75],[49.75,5.75],[50,5.75],[50.25,5.75],[50.5,5.75],[50.75,5.75],[51,5.75],[51.25,5.75],[51.5,5.75],[51.75,5.75],[52,5.75],[52.25,5.75],[52.5,5.75],[52.75,5.75],[53,5.75],[53.25,5.75],[45.5,6],[45.75,6],[46,6],[46.25,6],[46.5,6],[46.75,6],[47,6],[47.25,6],[47.5,6],[47.75,6],[48,6],[48.25,6],[48.5,6],[48.75,6],[49,6],[49.25,6],[49.5,6],[49.75,6],[50,6],[50.25,6],[50.5,6],[50.75,6],[51,6],[51.25,6],[51.5,6],[51.75,6],[52,6],[52.25,6],[52.5,6],[52.75,6],[53,6],[53.25,6],[45.75,6.25],[46,6.25],[46.25,6.25],[46.5,6.25],[46.75,6.25],[47,6.25],[47.25,6.25],[47.5,6.25],[47.75,6.25],[48,6.25],[48.25,6.25],[48.5,6.25],[48.75,6.25],[49,6.25],[49.25,6.25],[49.5,6.25],[49.75,6.25],[50,6.25],[50.25,6.25],[50.5,6.25],[50.75,6.25],[51,6.25],[51.25,6.25],[51.5,6.25],[51.75,6.25],[52,6.25],[52.25,6.25],[52.5,6.25],[52.75,6.25],[53,6.25],[53.25,6.25],[46,6.5],[46.25,6.5],[46.5,6.5],[46.75,6.5],[47,6.5],[47.25,6.5],[47.5,6.5],[47.75,6.5],[48,6.5],[48.25,6.5],[48.5,6.5],[48.75,6.5],[49,6.5],[49.25,6.5],[49.5,6.5],[49.75,6.5],[50,6.5],[50.25,6.5],[50.5,6.5],[50.75,6.5],[51,6.5],[51.25,6.5],[51.5,6.5],[51.75,6.5],[52,6.5],[52.25,6.5],[52.5,6.5],[52.75,6.5],[53,6.5],[53.25,6.5],[46.25,6.75],[46.5,6.75],[46.75,6.75],[47,6.75],[47.25,6.75],[47.5,6.75],[47.75,6.75],[48,6.75],[48.25,6.75],[48.5,6.75],[48.75,6.75],[49,6.75],[49.25,6.75],[49.5,6.75],[49.75,6.75],[50,6.75],[50.25,6.75],[50.5,6.75],[50.75,6.75],[51,6.75],[51.25,6.75],[51.5,6.75],[51.75,6.75],[52,6.75],[52.25,6.75],[52.5,6.75],[52.75,6.75],[53,6.75],[53.25,6.75],[46.25,7],[46.5,7],[46.75,7],[47,7],[47.25,7],[47.5,7],[47.75,7],[48,7],[48.25,7],[48.5,7],[48.75,7],[49,7],[49.25,7],[49.5,7],[49.75,7],[50,7],[50.25,7],[50.5,7],[50.75,7],[51,7],[51.25,7],[51.5,7],[51.75,7],[52,7],[52.25,7],[52.5,7],[52.75,7],[53,7],[53.25,7],[46.25,7.25],[46.5,7.25],[46.75,7.25],[47,7.25],[47.25,7.25],[47.5,7.25],[47.75,7.25],[48,7.25],[48.25,7.25],[48.5,7.25],[48.75,7.25],[49,7.25],[49.25,7.25],[49.5,7.25],[49.75,7.25],[50,7.25],[50.25,7.25],[50.5,7.25],[50.75,7.25],[51,7.25],[51.25,7.25],[51.5,7.25],[51.75,7.25],[52,7.25],[52.25,7.25],[52.5,7.25],[52.75,7.25],[53,7.25],[53.25,7.25],[53.5,7.25],[46.75,7.5],[47,7.5],[47.25,7.5],[47.5,7.5],[47.75,7.5],[48,7.5],[48.25,7.5],[48.5,7.5],[48.75,7.5],[49,7.5],[49.25,7.5],[49.5,7.5],[49.75,7.5],[50,7.5],[50.25,7.5],[50.5,7.5],[50.75,7.5],[51,7.5],[51.25,7.5],[51.5,7.5],[51.75,7.5],[52,7.5],[52.25,7.5],[52.5,7.5],[52.75,7.5],[53,7.5],[53.25,7.5],[53.5,7.5],[46.75,7.75],[47,7.75],[47.25,7.75],[47.5,7.75],[47.75,7.75],[48,7.75],[48.25,7.75],[48.5,7.75],[48.75,7.75],[49,7.75],[49.25,7.75],[49.5,7.75],[49.75,7.75],[50,7.75],[50.25,7.75],[50.5,7.75],[50.75,7.75],[51,7.75],[51.25,7.75],[51.5,7.75],[51.75,7.75],[52,7.75],[52.25,7.75],[52.5,7.75],[52.75,7.75],[53,7.75],[53.25,7.75],[53.5,7.75],[46.75,8],[47,8],[47.25,8],[47.5,8],[47.75,8],[48,8],[48.25,8],[48.5,8],[48.75,8],[49,8],[49.25,8],[49.5,8],[49.75,8],[50,8],[50.25,8],[50.5,8],[50.75,8],[51,8],[51.25,8],[51.5,8],[51.75,8],[52,8],[52.25,8],[52.5,8],[52.75,8],[53,8],[53.25,8],[53.5,8],[46.75,8.25],[47,8.25],[47.25,8.25],[47.5,8.25],[47.75,8.25],[48,8.25],[48.25,8.25],[48.5,8.25],[48.75,8.25],[49,8.25],[49.25,8.25],[49.5,8.25],[49.75,8.25],[50,8.25],[50.25,8.25],[50.5,8.25],[50.75,8.25],[51,8.25],[51.25,8.25],[51.5,8.25],[51.75,8.25],[52,8.25],[52.25,8.25],[52.5,8.25],[52.75,8.25],[53,8.25],[53.25,8.25],[53.5,8.25],[46.75,8.5],[47,8.5],[47.25,8.5],[47.5,8.5],[47.75,8.5],[48,8.5],[48.25,8.5],[48.5,8.5],[48.75,8.5],[49,8.5],[49.25,8.5],[49.5,8.5],[49.75,8.5],[50,8.5],[50.25,8.5],[50.5,8.5],[50.75,8.5],[51,8.5],[51.25,8.5],[51.5,8.5],[51.75,8.5],[52,8.5],[52.25,8.5],[52.5,8.5],[52.75,8.5],[53,8.5],[53.25,8.5],[53.5,8.5],[47,8.75],[47.25,8.75],[47.5,8.75],[47.75,8.75],[48,8.75],[48.25,8.75],[48.5,8.75],[48.75,8.75],[49,8.75],[49.25,8.75],[49.5,8.75],[49.75,8.75],[50,8.75],[50.25,8.75],[50.5,8.75],[50.75,8.75],[51,8.75],[51.25,8.75],[51.5,8.75],[51.75,8.75],[52,8.75],[52.25,8.75],[52.5,8.75],[52.75,8.75],[53,8.75],[53.25,8.75],[53.5,8.75],[53.75,8.75],[54.75,8.75],[55,8.75],[47.25,9],[47.5,9],[47.75,9],[48,9],[48.25,9],[48.5,9],[48.75,9],[49,9],[49.25,9],[49.5,9],[49.75,9],[50,9],[50.25,9],[50.5,9],[50.75,9],[51,9],[51.25,9],[51.5,9],[51.75,9],[52,9],[52.25,9],[52.5,9],[52.75,9],[53,9],[53.25,9],[53.5,9],[53.75,9],[54,9],[54.25,9],[54.5,9],[54.75,9],[55,9],[47.25,9.25],[47.5,9.25],[47.75,9.25],[48,9.25],[48.25,9.25],[48.5,9.25],[48.75,9.25],[49,9.25],[49.25,9.25],[49.5,9.25],[49.75,9.25],[50,9.25],[50.25,9.25],[50.5,9.25],[50.75,9.25],[51,9.25],[51.25,9.25],[51.5,9.25],[51.75,9.25],[52,9.25],[52.25,9.25],[52.5,9.25],[52.75,9.25],[53,9.25],[53.25,9.25],[53.5,9.25],[53.75,9.25],[54,9.25],[54.25,9.25],[54.5,9.25],[54.75,9.25],[55,9.25],[47.25,9.5],[47.5,9.5],[47.75,9.5],[48,9.5],[48.25,9.5],[48.5,9.5],[48.75,9.5],[49,9.5],[49.25,9.5],[49.5,9.5],[49.75,9.5],[50,9.5],[50.25,9.5],[50.5,9.5],[50.75,9.5],[51,9.5],[51.25,9.5],[51.5,9.5],[51.75,9.5],[52,9.5],[52.25,9.5],[52.5,9.5],[52.75,9.5],[53,9.5],[53.25,9.5],[53.5,9.5],[53.75,9.5],[54,9.5],[54.25,9.5],[54.5,9.5],[54.75,9.5],[47.25,9.75],[47.5,9.75],[47.75,9.75],[48,9.75],[48.25,9.75],[48.5,9.75],[48.75,9.75],[49,9.75],[49.25,9.75],[49.5,9.75],[49.75,9.75],[50,9.75],[50.25,9.75],[50.5,9.75],[50.75,9.75],[51,9.75],[51.25,9.75],[51.5,9.75],[51.75,9.75],[52,9.75],[52.25,9.75],[52.5,9.75],[52.75,9.75],[53,9.75],[53.25,9.75],[53.5,9.75],[53.75,9.75],[54,9.75],[54.25,9.75],[54.5,9.75],[54.75,9.75],[47.25,10],[47.5,10],[47.75,10],[48,10],[48.25,10],[48.5,10],[48.75,10],[49,10],[49.25,10],[49.5,10],[49.75,10],[50,10],[50.25,10],[50.5,10],[50.75,10],[51,10],[51.25,10],[51.5,10],[51.75,10],[52,10],[52.25,10],[52.5,10],[52.75,10],[53,10],[53.25,10],[53.5,10],[53.75,10],[54,10],[54.25,10],[54.5,10],[54.75,10],[47.25,10.25],[47.5,10.25],[47.75,10.25],[48,10.25],[48.25,10.25],[48.5,10.25],[48.75,10.25],[49,10.25],[49.25,10.25],[49.5,10.25],[49.75,10.25],[50,10.25],[50.25,10.25],[50.5,10.25],[50.75,10.25],[51,10.25],[51.25,10.25],[51.5,10.25],[51.75,10.25],[52,10.25],[52.25,10.25],[52.5,10.25],[52.75,10.25],[53,10.25],[53.25,10.25],[53.5,10.25],[53.75,10.25],[54,10.25],[54.25,10.25],[47.5,10.5],[47.75,10.5],[48,10.5],[48.25,10.5],[48.5,10.5],[48.75,10.5],[49,10.5],[49.25,10.5],[49.5,10.5],[49.75,10.5],[50,10.5],[50.25,10.5],[50.5,10.5],[50.75,10.5],[51,10.5],[51.25,10.5],[51.5,10.5],[51.75,10.5],[52,10.5],[52.25,10.5],[52.5,10.5],[52.75,10.5],[53,10.5],[53.25,10.5],[53.5,10.5],[53.75,10.5],[54,10.5],[54.25,10.5],[47.5,10.75],[47.75,10.75],[48,10.75],[48.25,10.75],[48.5,10.75],[48.75,10.75],[49,10.75],[49.25,10.75],[49.5,10.75],[49.75,10.75],[50,10.75],[50.25,10.75],[50.5,10.75],[50.75,10.75],[51,10.75],[51.25,10.75],[51.5,10.75],[51.75,10.75],[52,10.75],[52.25,10.75],[52.5,10.75],[52.75,10.75],[53,10.75],[53.25,10.75],[53.5,10.75],[53.75,10.75],[54,10.75],[54.25,10.75],[47.5,11],[47.75,11],[48,11],[48.25,11],[48.5,11],[48.75,11],[49,11],[49.25,11],[49.5,11],[49.75,11],[50,11],[50.25,11],[50.5,11],[50.75,11],[51,11],[51.25,11],[51.5,11],[51.75,11],[52,11],[52.25,11],[52.5,11],[52.75,11],[53,11],[53.25,11],[53.5,11],[53.75,11],[54,11],[47.5,11.25],[47.75,11.25],[48,11.25],[48.25,11.25],[48.5,11.25],[48.75,11.25],[49,11.25],[49.25,11.25],[49.5,11.25],[49.75,11.25],[50,11.25],[50.25,11.25],[50.5,11.25],[50.75,11.25],[51,11.25],[51.25,11.25],[51.5,11.25],[51.75,11.25],[52,11.25],[52.25,11.25],[52.5,11.25],[52.75,11.25],[53,11.25],[53.25,11.25],[53.5,11.25],[53.75,11.25],[54,11.25],[47.5,11.5],[47.75,11.5],[48,11.5],[48.25,11.5],[48.5,11.5],[48.75,11.5],[49,11.5],[49.25,11.5],[49.5,11.5],[49.75,11.5],[50,11.5],[50.25,11.5],[50.5,11.5],[50.75,11.5],[51,11.5],[51.25,11.5],[51.5,11.5],[51.75,11.5],[52,11.5],[52.25,11.5],[52.5,11.5],[52.75,11.5],[53,11.5],[53.25,11.5],[53.5,11.5],[53.75,11.5],[54,11.5],[47.5,11.75],[47.75,11.75],[48,11.75],[48.25,11.75],[48.5,11.75],[48.75,11.75],[49,11.75],[49.25,11.75],[49.5,11.75],[49.75,11.75],[50,11.75],[50.25,11.75],[50.5,11.75],[50.75,11.75],[51,11.75],[51.25,11.75],[51.5,11.75],[51.75,11.75],[52,11.75],[52.25,11.75],[52.5,11.75],[52.75,11.75],[53,11.75],[53.25,11.75],[53.5,11.75],[53.75,11.75],[54,11.75],[47.5,12],[47.75,12],[48,12],[48.25,12],[48.5,12],[48.75,12],[49,12],[49.25,12],[49.5,12],[49.75,12],[50,12],[50.25,12],[50.5,12],[50.75,12],[51,12],[51.25,12],[51.5,12],[51.75,12],[52,12],[52.25,12],[52.5,12],[52.75,12],[53,12],[53.25,12],[53.5,12],[53.75,12],[54,12],[47.5,12.25],[47.75,12.25],[48,12.25],[48.25,12.25],[48.5,12.25],[48.75,12.25],[49,12.25],[49.25,12.25],[49.5,12.25],[49.75,12.25],[50,12.25],[50.25,12.25],[50.5,12.25],[50.75,12.25],[51,12.25],[51.25,12.25],[51.5,12.25],[51.75,12.25],[52,12.25],[52.25,12.25],[52.5,12.25],[52.75,12.25],[53,12.25],[53.25,12.25],[53.5,12.25],[53.75,12.25],[54,12.25],[54.25,12.25],[47.5,12.5],[47.75,12.5],[48,12.5],[48.25,12.5],[48.5,12.5],[48.75,12.5],[49,12.5],[49.25,12.5],[49.5,12.5],[49.75,12.5],[50,12.5],[50.25,12.5],[50.5,12.5],[50.75,12.5],[51,12.5],[51.25,12.5],[51.5,12.5],[51.75,12.5],[52,12.5],[52.25,12.5],[52.5,12.5],[52.75,12.5],[53,12.5],[53.25,12.5],[53.5,12.5],[53.75,12.5],[54,12.5],[54.25,12.5],[47.5,12.75],[47.75,12.75],[48,12.75],[48.25,12.75],[48.5,12.75],[48.75,12.75],[49,12.75],[49.25,12.75],[49.5,12.75],[49.75,12.75],[50,12.75],[50.25,12.75],[50.5,12.75],[50.75,12.75],[51,12.75],[51.25,12.75],[51.5,12.75],[51.75,12.75],[52,12.75],[52.25,12.75],[52.5,12.75],[52.75,12.75],[53,12.75],[53.25,12.75],[53.5,12.75],[53.75,12.75],[54,12.75],[54.25,12.75],[47.5,13],[47.75,13],[48,13],[48.25,13],[48.5,13],[48.75,13],[49,13],[49.25,13],[49.5,13],[49.75,13],[50,13],[50.25,13],[50.5,13],[50.75,13],[51,13],[51.25,13],[51.5,13],[51.75,13],[52,13],[52.25,13],[52.5,13],[52.75,13],[53,13],[53.25,13],[53.5,13],[53.75,13],[54,13],[54.25,13],[47.5,13.25],[47.75,13.25],[48,13.25],[48.25,13.25],[48.5,13.25],[48.75,13.25],[49,13.25],[49.25,13.25],[49.5,13.25],[49.75,13.25],[50,13.25],[50.25,13.25],[50.5,13.25],[50.75,13.25],[51,13.25],[51.25,13.25],[51.5,13.25],[51.75,13.25],[52,13.25],[52.25,13.25],[52.5,13.25],[52.75,13.25],[53,13.25],[53.25,13.25],[53.5,13.25],[53.75,13.25],[54,13.25],[54.25,13.25],[54.5,13.25],[47.5,13.5],[47.75,13.5],[48,13.5],[48.25,13.5],[48.5,13.5],[48.75,13.5],[49,13.5],[49.25,13.5],[49.5,13.5],[49.75,13.5],[50,13.5],[50.25,13.5],[50.5,13.5],[50.75,13.5],[51,13.5],[51.25,13.5],[51.5,13.5],[51.75,13.5],[52,13.5],[52.25,13.5],[52.5,13.5],[52.75,13.5],[53,13.5],[53.25,13.5],[53.5,13.5],[53.75,13.5],[54,13.5],[54.25,13.5],[54.5,13.5],[47.5,13.75],[47.75,13.75],[48,13.75],[48.25,13.75],[48.5,13.75],[48.75,13.75],[49,13.75],[49.25,13.75],[49.5,13.75],[49.75,13.75],[50,13.75],[50.25,13.75],[50.5,13.75],[50.75,13.75],[51,13.75],[51.25,13.75],[51.5,13.75],[51.75,13.75],[52,13.75],[52.25,13.75],[52.5,13.75],[52.75,13.75],[53,13.75],[53.25,13.75],[53.5,13.75],[53.75,13.75],[54,13.75],[47.75,14],[48,14],[48.25,14],[48.5,14],[48.75,14],[49,14],[49.25,14],[49.5,14],[49.75,14],[50,14],[50.25,14],[50.5,14],[50.75,14],[51,14],[51.25,14],[51.5,14],[51.75,14],[52,14],[52.25,14],[52.5,14],[52.75,14],[53,14],[53.25,14],[53.5,14],[53.75,14],[49,14.25],[49.25,14.25],[49.5,14.25],[49.75,14.25],[50,14.25],[50.25,14.25],[50.5,14.25],[50.75,14.25],[51,14.25],[51.25,14.25],[51.5,14.25],[51.75,14.25],[52,14.25],[52.25,14.25],[52.5,14.25],[52.75,14.25],[53,14.25],[53.25,14.25],[53.5,14.25],[53.75,14.25],[50,14.5],[50.25,14.5],[50.5,14.5],[50.75,14.5],[51,14.5],[51.25,14.5],[51.5,14.5],[51.75,14.5],[52,14.5],[52.25,14.5],[52.5,14.5],[52.75,14.5],[53,14.5],[53.25,14.5],[53.5,14.5],[53.75,14.5],[50,14.75],[50.25,14.75],[50.5,14.75],[50.75,14.75],[51,14.75],[51.25,14.75],[51.5,14.75],[51.75,14.75],[52,14.75],[52.25,14.75],[52.5,14.75],[52.75,14.75],[53,14.75],[53.25,14.75],[53.5,14.75],[53.75,14.75],[50.25,15],[50.5,15],[50.75,15],[51,15],[51.25,15],[51.5,15],[51.75,15],[52,15],[52.25,15],[52.5,15],[52.75,15],[53,15],[53.25,15],[53.5,15],[53.75,15],[50.25,15.25],[50.5,15.25],[50.75,15.25],[51,15.25],[51.25,15.25],[51.5,15.25],[51.75,15.25],[52,15.25],[52.25,15.25],[52.5,15.25],[52.75,15.25],[53,15.25],[53.25,15.25],[53.5,15.25],[50.5,15.5],[50.75,15.5],[51,15.5],[51.25,15.5],[51.5,15.5],[51.75,15.5],[52,15.5],[52.25,15.5],[52.5,15.5],[52.75,15.5],[53,15.5],[53.25,15.5],[53.5,15.5],[50.75,15.75],[51,15.75],[51.25,15.75],[51.5,15.75],[51.75,15.75],[52,15.75],[52.25,15.75],[52.5,15.75],[52.75,15.75],[53,15.75],[53.25,15.75],[52.5,16],[52.75,16],[53,16]]

t0 = datetime.date(2018, 1, 1)

d=[]
for ll in latlon:
    d2=[]
    a = astral.Location(('name', 'region', ll[0], ll[1], 'UTC', 0))
    for t in range(0,366):
        try:
            tmp =a.sun(t0 + datetime.timedelta(t))
            d2.append([str(tmp['dawn']).split('+')[0], str(tmp['sunrise']).split('+')[0], str(tmp['sunset']).split('+')[0], str(tmp['dusk']).split('+')[0]])
        except:
            d2.append([])
    d.append(d2)


with open('sunriseset_grid.json', 'w') as fout:
    json.dump(d, fout)
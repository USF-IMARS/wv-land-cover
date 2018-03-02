import os, string, sys, shutil, glob, re, tarfile, logging, argparse, signal
from datetime import *
from subprocess import *
from math import *
from xml.etree import cElementTree as ET

import gdal, ogr,osr, gdalconst
import numpy

logger = logging.getLogger("logger")
logger.setLevel(logging.DEBUG)


MODES = ["ALL","MOSAIC","SHP","TEST"]
EXTS = [".tif"]
GTIFF_COMPRESSIONS = ["jpeg95","lzw"]

#class Attribs:
#    def __init__(self,dAttribs):
#        self.cc = dAttribs["cc"]
#        self.sunel = dAttribs["sunel"]
#        self.ona = dAttribs["ona"]
#        self.tdi = dAttribs["tdi"]
#        self.alr = dAttribs["alr"]
#        self.exdur = dAttribs["exdur"]
#        self.datediff = dAttribs["datediff"]
#        self.exfact = dAttribs["exfact"]
#        self.panfact = dAttribs["panfact"]


def build_arg_list(args,pos_arg_keys,arg_keys_to_remove):
    args_dict = vars(args)
    arg_list = []
    
    ## Add optional args to arg_list
    for k,v in args_dict.iteritems():
        if k not in pos_arg_keys and k not in arg_keys_to_remove and v is not None:
            if isinstance(v,list) or isinstance(v,tuple):
                arg_list.append("--%s %s" %(k,' '.join([str(item) for item in v])))
            elif isinstance(v,bool):
                if v is True:
                    arg_list.append("--%s" %(k))
            else:
                arg_list.append("--%s %s" %(k,str(v)))
    
    arg_str = " ".join(arg_list)
    return arg_str


def buildMosaicParentArgumentParser():
    
    #### Set Up Arguments 
    parser = argparse.ArgumentParser(add_help=False)
    
    ####Optional Arguments
    
    parser.add_argument("-r", "--resolution", nargs=2, type=float,
                        help="output pixel resolution -- xres yres (default is same as first input file)")
    parser.add_argument("-e", "--extent", nargs=4, type=float,
                        help="extent of output mosaic -- xmin xmax ymin ymax (default is union of all inputs)")
    parser.add_argument("-t", "--tilesize", nargs=2, type=float,
                        help="tile size in coordinate system units -- xsize ysize (default is 40,000 times output resolution)")
    parser.add_argument("--force_pan_to_multi", action="store_true", dest="force_pan_to_multi", default=False,
                        help="if output is multiband, force script to also use 1 band images")
    parser.add_argument("-b", "--bands", type=int,
                        help="number of output bands( default is number of bands in the first image)")
    parser.add_argument("--tday",
                        help="month and day of the year to use as target for image suitability ranking -- 04-05")
    parser.add_argument("--nosort", action="store_true", default=False,
                        help="do not sort images by metadata. script uses the order of the input textfile or directory (first image is first drawn).  Not recommended if input is a directory; order will be random")
    parser.add_argument("--use_exposure", action="store_true", default=False,
                        help="use exposure settings in metadata to inform score")
    parser.add_argument("--exclude",
                        help="file of file name patterns (text only, no wildcards or regexs) to exclude")

    return parser


class ImageInfo:
    def __init__(self,src,frmt,srs=None):
        
        self.frmt = frmt  #image format (IMAGE,RECORD)
        
        if frmt == 'IMAGE':
            self.get_attributes_from_file(src)
        elif frmt == 'RECORD':
            self.get_attributes_from_record(src,srs)
        else:
            logger.error("Image format must be RECORD or IMAGE")
        
        
    #self.xsize = None
    #self.ysize = None
    #self.proj = None
    #self.bands = None
    #self.datatype = None
    #self.datatype_readable = None
    #self.xres = None
    #self.yres = None
    #self.geom = None
    #self.sensor = None
    #self.acqdate = None
    #"cc":None,
    #"sunel":None,
    #"ona":None,
    #"date":None,
    #"tdi":None
        
    def get_attributes_from_record(self, feat, srs):
                
        i = feat.GetFieldIndex("S_FILEPATH")
        if i == -1:
            i = feat.GetFieldIndex("O_FILEPATH")
        elif len(feat.GetFieldAsString(i)) == 0:
            i = feat.GetFieldIndex("O_FILEPATH")
        path = feat.GetFieldAsString(i)
        
        if len(path) > 1:
            if r"V:/pgc/agic/private" in path:
                srcfp = path.replace(r"V:/pgc",r'/mnt/agic/storage00')
            elif r"/pgc/agic/private" in path:
                srcfp = path.replace(r"/pgc",r'/mnt/agic/storage00')
            elif r"V:/pgc/data" in path:
                srcfp = path.replace(r"V:/pgc/data",r'/mnt/pgc/data')
            elif r"/pgc/data" in path:
                srcfp = path.replace(r"/pgc/data",r'/mnt/pgc/data')
            else:
                srcfp = path
            
        self.srcfp = srcfp
        self.srcdir, self.srcfn = os.path.split(srcfp)
        
        i = feat.GetFieldIndex("COLUMNS")
        if i != -1:
            self.xsize = feat.GetFieldAsDouble(i)
        i = feat.GetFieldIndex("ROWS")
        if i != -1:
            self.ysize = feat.GetFieldAsDouble(i)
        i = feat.GetFieldIndex("BANDS")
        if i != -1:
            self.bands = feat.GetFieldAsDouble(i)
        
        self.proj = srs.ExportToWkt()
        self.xres = None
        self.yres = None
        self.datatype = None
        
        i = feat.GetFieldIndex("SUN_ELEV")
        if i != -1:
            self.sunel = feat.GetFieldAsDouble(i)
        i = feat.GetFieldIndex("OFF_NADIR")
        if i != -1:
            self.ona = feat.GetFieldAsDouble(i)
        i = feat.GetFieldIndex("CLOUDCOVER")
        if i != -1:
            self.cloudcover = feat.GetFieldAsDouble(i)
        i = feat.GetFieldIndex("SENSOR")
        if i != -1:
            self.sensor = feat.GetFieldAsString(i)
        i = feat.GetFieldIndex("SCENE_ID")
        if i != -1:
            self.scene_id = feat.GetFieldAsString(i)
        i = feat.GetFieldIndex("CATALOG_ID")
        if i != -1:
            self.catid = feat.GetFieldAsString(i)
        
        i = feat.GetFieldIndex("TDI")
        if i != -1:
            tdi_str = feat.GetFieldAsString(i)
            tdi_list = tdi_str.split('|')
            self.tdi = None
            for item in tdi_list:
                if 'pan' in item:
                    self.tdi = int(item[4:])
                if 'green' in item:
                    self.tdi = int(item[6:])
        
        i = feat.GetFieldIndex("ACQ_TIME")
        if i != -1:
            date_str = feat.GetFieldAsString(i)
            self.acqdate = datetime.strptime(date_str[:19],"%Y-%m-%dT%H:%M:%S")
        
        geom = feat.GetGeometryRef()
        self.geom = geom.Clone()
        
    
    def get_attributes_from_file(self, srcfp):
        self.srcfp = srcfp
        self.srcdir, self.srcfn = os.path.split(srcfp)
        
        ds = gdal.Open(self.srcfp)
        if ds is not None:
            self.xsize = ds.RasterXSize
            self.ysize = ds.RasterYSize
            self.proj = ds.GetProjectionRef() if ds.GetProjectionRef() != '' else ds.GetGCPProjection()
            self.bands = ds.RasterCount
            self.datatype = ds.GetRasterBand(1).DataType
            self.datatype_readable = gdal.GetDataTypeName(self.datatype)

            gtf = ds.GetGeoTransform()
            num_gcps = ds.GetGCPCount()
            
            if num_gcps == 0:
                
                self.xres = abs(gtf[1])
                self.yres = abs(gtf[5])
                ulx = gtf[0] + 0 * gtf[1] + 0 * gtf[2]
                uly = gtf[3] + 0 * gtf[4] + 0 * gtf[5]
                urx = gtf[0] + self.xsize * gtf[1] + 0 * gtf[2]
                ury = gtf[3] + self.xsize * gtf[4] + 0 * gtf[5]
                llx = gtf[0] + 0 * gtf[1] + self.ysize * gtf[2]
                lly = gtf[3] + 0 * gtf[4] + self.ysize * gtf[5]
                lrx = gtf[0] + self.xsize * gtf[1] + self.ysize* gtf[2]
                lry = gtf[3] + self.xsize * gtf[4] + self.ysize * gtf[5]
                
            
            elif num_gcps == 4:
                
                gcps = ds.GetGCPs()
                gcp_dict = {}
                id_dict = {"UpperLeft":1,
                           "1":1,
                           "UpperRight":2,
                           "2":2,
                           "LowerLeft":4,
                           "4":4,
                           "LowerRight":3,
                           "3":3}
        
                for gcp in gcps:
                    gcp_dict[id_dict[gcp.Id]] = [float(gcp.GCPPixel), float(gcp.GCPLine), float(gcp.GCPX), float(gcp.GCPY), float(gcp.GCPZ)]
        
                ulx = gcp_dict[1][2]
                uly = gcp_dict[1][3]
                urx = gcp_dict[2][2]
                ury = gcp_dict[2][3]
                llx = gcp_dict[4][2]
                lly = gcp_dict[4][3]
                lrx = gcp_dict[3][2]
                lry = gcp_dict[3][3]
        
                self.xres = abs(math.sqrt((ulx - urx)**2 + (uly - ury)**2)/ self.xsize)
                self.yres = abs(math.sqrt((ulx - llx)**2 + (uly - lly)**2)/ self.ysize)
                
            poly_wkt = 'POLYGON (( %.12f %.12f, %.12f %.12f, %.12f %.12f, %.12f %.12f, %.12f %.12f ))' %(ulx,uly,urx,ury,lrx,lry,llx,lly,ulx,uly)
            self.geom = ogr.CreateGeometryFromWkt(poly_wkt)
            self.xs = [ulx,urx,lrx,llx]
            self.ys = [uly,ury,lry,lly]
                
                
        else:
            logger.warning("Cannot open image: %s" %self.srcfp)
            self.xsize = None
            self.ysize = None
            self.proj = None
            self.bands = None
            self.datatype = None
            self.datatype_readable = None
            self.xres = None
            self.yres = None

        ds = None
        
        #### Set unknown attribs to None for now
        self.sunel = None
        self.ona = None
        self.cloudcover = None
        self.sensor = None
        self.scene_id = None
        self.catid = None
        self.tdi = None
        self.acqdate = None
 
   
    def get_attributes_from_xml(self):
 
        
        dAttribs = {
            "cc":None,
            "sunel":None,
            "ona":None,
            "date":None,
            "tdi":None,
            "catid":None,
            "sensor":None
        }
    
        dTags = {
            ## DG tags
            "CATID":"catid",
            "SATID":"sensor",
            "CLOUDCOVER":"cc",
            "MEANSUNEL":"sunel",
            "MEANOFFNADIRVIEWANGLE":"ona",
            "FIRSTLINETIME":"date",
            "TDILEVEL":"tdi",
            
            ## GE tags
            "archiveId":"catid",
            "satelliteName":"sensor",
            "percentCloudCover":"cc",
            "firstLineSunElevationAngle":"sunel",
            "firstLineElevationAngle":"ona",
            "firstLineAcquisitionDateTime":"date",
            "tdiMode":"tdi"
        }
            
        paths = (
            os.path.splitext(self.srcfp)[0]+'.xml',
            os.path.splitext(self.srcfp)[0]+'.XML',
            os.path.splitext(self.srcfp)[0]+'.txt',
            os.path.splitext(self.srcfp)[0]+'.pvl',
        )
        
        metapath = None
        for path in paths:
            if os.path.isfile(path):
                metapath = path
                break
        
        if not metapath:
            logger.debug("No metadata found for %s" %self.srcfp)
        
        else:
            metad = None
            
            #### if xml format
            if os.path.splitext(metapath)[1].lower() == '.xml':
                try:
                    metad = ET.parse(metapath)
                except ET.ParseError, err:
                    logger.debug("ERROR parsing metadata: %s, %s" %(err,metapath))
            
            else:
                try:
                    metad = getGEMetadataAsXml(metapath)
                except Exception, err:
                    logger.debug("ERROR parsing metadata: %s, %s" %(err,metapath))
                #### Write IK01 code 
        
            if metad is not None:
                
                for tag in dTags:
                    taglist = metad.findall(".//%s"%tag)
                    vallist = []
                    for elem in taglist:
                        
                        text = elem.text
                    
                        if text is not None:
                            try:
                                if tag == "firstLineElevationAngle":
                                    val = 90 - float(text)
                                elif tag in ["FIRSTLINETIME","firstLineAcquisitionDateTime","CATID","archiveId","SATID"]:
                                    val = text
                                elif tag == "percentCloudCover":
                                    val = float(text)/100
                                elif tag == "satelliteName":
                                    val = "GE01"
                                else:
                                    val = float(text)
                                    
                                vallist.append(val)
                                
                            except Exception, e:
                                logger.debug("Error reading metadata values: %s, %s" %(metapath,e))
                                
                    if dTags[tag] == 'tdi' and len(taglist) > 1:    
                        #### use pan or green band TDI for exposure calculation
                        if len(vallist) == 4:
                            dAttribs['tdi'] = vallist[1]
                        elif len(vallist) == 5 and self.bands == 1: #pan image
                            dAttribs['tdi'] = vallist[4]
                        elif len(vallist) == 5 and self.bands in [3,4]: #multi image
                            dAttribs['tdi'] = vallist[1]
                        elif len(vallist) == 8:
                            dAttribs['tdi'] = vallist[3]
                        else:
                            logger.debug("Unexpected number of TDI values and band count ( TDI: expected 1, 4, 5, or 8 - found %d ; Band cound, expected 1, 4, or 8 - found %d) %s" %(len(vallist), self.bands, metapath))
                            
                    elif len(taglist) == 1:
                        val = vallist[0]
                        dAttribs[dTags[tag]] = val
                        
                    elif len(taglist) <> 0:
                        logger.debug("Unexpected number of %s values, %s" %(tag,metapath))
                
                self.sunel = dAttribs["sunel"]
                self.ona = dAttribs["ona"]
                self.cloudcover = dAttribs["cc"]
                self.sensor = dAttribs["sensor"]
                self.catid = dAttribs["catid"]
                self.tdi = dAttribs["tdi"]
                self.acqdate = datetime.strptime(dAttribs["date"],"%Y-%m-%dT%H:%M:%S.%fZ")
                

    def getScore(self,params):
        
        score = 0
       
        if not self.catid:
            self.get_attributes_from_xml()
        
        required_attribs = [
            self.sunel,
            self.ona,
            self.cloudcover,
            self.sensor,
        ]
        
        #### Test if all required values were found in metadata search
        status = [val is None for val in required_attribs]
        
        if sum(status) != 0:
            logger.error("Cannot determine score for image {0}:\n  Sun elev\t{1}\n  Off nadir\t{2}\n  Cloudcover\t{3}\n  Sensor\t{4}".format(self.srcfn,self.sunel,self.ona,self.cloudcover,self.sensor))
            score = -1
        
        else:
            
            #### Assign panfactor if pan images are to be included in a multispectral mosaic   
            if self.bands == 1 and params.force_pan_to_multi is True:
                self.panfactor = 0.5
            else:
                self.panfactor = 1
                    
            #### Test if TDI is needed, get exposure factor
            if params.useExposure is True:
                if self.tdi is None:
                    logger.error("Cannot get tdi for image to determine exposure settings: {0}".format(self.srcfn))
                    self.exposure_factor = None
                else:
                    exfact = self.tdi * self.sunel
                    self.exposure_factor = exfact
                    
                    pan_exposure_thresholds = {
                        "WV01":1400,
                        "WV02":1400,
                        "WV03":1400,
                        "QB02":500,
                        #"GE01":,
                    }
                    
                    multi_exposure_thresholds = {
                        "WV02":400,
                        "WV03":400,
                        "GE01":170,
                        "QB02":25,
                    }
                    
                    #### Remove images with high exposure settings (tdi_pan (or tdi_grn) * sunel)
                    if params.bands == 1:
                        if self.sensor in pan_exposure_thresholds:
                            if exfact > pan_exposure_thresholds[self.sensor]:
                                logger.debug("Image overexposed: %s --> %i" %(self.srcfp,exfact))
                                score = -1
                    
                    else:
                        if self.sensor in multi_exposure_thresholds:
                            if exfact > multi_exposure_thresholds[self.sensor]:
                                logger.debug("Image overexposed: %s --> %i" %(self.srcfp,exfact))
                                score = -1
            
            #### Test if acqdate if needed, get date difference
            if params.m != 0:
                if self.acqdate is None:
                    logger.error("Cannot get acqdate for image to determine date-based score: {0}".format(self.srcfn))
                    self.date_diff = -9999
                    
                else:
                    #### Find nearest year for target day
                    tdeltas = []
                    for y in range(self.acqdate.year-1,self.acqdate.year+2):
                        tdeltas.append(abs((datetime(y,params.m,params.d) - self.acqdate).days))
                    
                    self.date_diff = min(tdeltas)
            
            
                #### Assign weights
                ccwt = 30
                sunelwt = 10
                onawt = 5
                datediffwt = 55
                
            else:
                self.date_diff = -9999
                ccwt = 48
                sunelwt = 28
                onawt = 24
                datediffwt = 0
                
                
            #### Handle nonesense or nodata cloud cover values
            if self.cloudcover < 0 or self.cloudcover > 1:
                self.cloudcover = 0.5
            
            if self.cloudcover > 0.5:
                logger.debug("Image too cloudy (>50 percent): %s --> %f" %(self.srcfp,self.cloudcover))
                score = -1
            
            #### Handle ridiculously low sun el values, these images will result is spurious TOA values
            if self.sunel < 2:
                logger.debug("Sun elevation too low (<2 degrees): %s --> %f" %(self.srcfp,self.sunel))
                score = -1
                        
            if not score == -1:
                rawscore = ccwt * (1-self.cloudcover) + sunelwt * (self.sunel/90) + onawt * ((90-self.ona)/90.0) + datediffwt * ((183 - self.date_diff)/183.0)
                score = rawscore * self.panfactor  
        
        self.score = score
        return self.score
    
        
class MosaicParams:
    pass

class TileParams:
    def __init__(self,x,x2,y,y2,j,i,name):
        self.xmin = x
        self.xmax = x2
        self.ymin = y
        self.ymax = y2
        self.i = i
        self.j = j
        self.name = name
        poly_wkt = 'POLYGON (( %f %f, %f %f, %f %f, %f %f, %f %f ))' %(x,y,x,y2,x2,y2,x2,y,x,y)
        self.geom = ogr.CreateGeometryFromWkt(poly_wkt)
        

def filterMatchingImages(imginfo_list,params):
    imginfo_list2 = []
    
    for iinfo in imginfo_list:
        #print iinfo.srcfp, iinfo.proj
        isSame = True
        p = osr.SpatialReference()
        p.ImportFromWkt(iinfo.proj)
        rp = osr.SpatialReference()
        rp.ImportFromWkt(params.proj)
        if p.IsSame(rp) is False:
            isSame = False
        if iinfo.bands != params.bands and not (params.force_pan_to_multi is True and iinfo.bands == 1):
            isSame = False
        if iinfo.datatype != params.datatype:
            isSame = False
            
        if isSame is True:
            imginfo_list2.append(iinfo)
        else:
            logger.debug("Image does not match filter: %s" %iinfo.srcfp)

    return imginfo_list2


def getMosaicParameters(iinfo,options):
    
    params = MosaicParams()
    
    if options.resolution is not None:
        params.xres = options.resolution[0]
        params.yres = options.resolution[1]
    else:
        params.xres = iinfo.xres
        params.yres = iinfo.yres
    
    params.bands = options.bands if options.bands is not None else iinfo.bands
    params.proj = iinfo.proj
    params.datatype = iinfo.datatype
    params.useExposure = options.use_exposure
    
    if options.tday is not None:
        params.m = int(options.tday.split("-")[0])
        params.d = int(options.tday.split("-")[1])   
    else:
        params.m = 0
        params.d = 0
    
    if options.extent is not None: # else set after geoms are collected
        params.xmin = options.extent[0]
        params.ymin = options.extent[2]
        params.xmax = options.extent[1]
        params.ymax = options.extent[3]
        
    if options.tilesize is not None:
        params.xtilesize = options.tilesize[0]
        params.ytilesize = options.tilesize[1]
    elif params.xres is not None:
        params.xtilesize = params.xres * 40000
        params.ytilesize = params.yres * 40000
    else:
        params.xtilesize = None
        params.ytilesize = None
    
    params.force_pan_to_multi = True if params.bands > 1 and options.force_pan_to_multi else False # determine if force pan to multi is applicable and true
    
    return params


def GetExactTrimmedGeom(image, step=2, tolerance=1):
    
    geom2 = None
    geom = None
    xs,ys = [],[]
    ds = gdal.Open(image)
    if ds is not None:
        if ds.RasterCount > 0:
            
            inband = ds.GetRasterBand(1)
            
            nd = inband.GetNoDataValue()
            if nd is None:
                nd = 0
            
            #print ("Image NoData Value: %d" %nd )
            
            gtf = ds.GetGeoTransform()
        
            pixelst = []
            pixelsb = []
            pts = []
            
            #### For every other line, find first and last data pixel
            lines = xrange(0, inband.YSize, step)
            
            xsize = inband.XSize
            npflatnonzero = numpy.flatnonzero
            bandReadAsArray = inband.ReadAsArray
            
            lines_flatnonzero = [npflatnonzero(bandReadAsArray(0,l,xsize,1) != nd) for l in lines]
            i = 0
            
            for nz in lines_flatnonzero:
                
                nzmin = nz[0] if nz.size > 0 else 0
                nzmax = nz[-1] if nz.size > 0 else 0
                
                if nz.size > 0:
                    pixelst.append((nzmax+1,i))
                    pixelsb.append((nzmin,i))           
                i += step
                
            pixelsb.reverse()
            pixels = pixelst + pixelsb
            
            #print len(pixels)
            
            for px in pixels:
                x,y = pl2xy(gtf,inband,px[0],px[1])
                xs.append(x)
                ys.append(y)
                pts.append((x,y))
                #print px[0],px[1],x,y
            
            #### create geometry
            poly_vts = []
            for pt in pts:
                poly_vts.append("%.16f %.16f" %(pt[0],pt[1]))
            if len(pts) > 0:
                poly_vts.append("%.16f %.16f" %(pts[0][0],pts[0][1]))
            
            if len(poly_vts) > 0:
                poly_wkt = 'POLYGON (( %s ))' %(string.join(poly_vts,", "))
                #print poly_wkt
                
                geom = ogr.CreateGeometryFromWkt(poly_wkt)
                #print geom
                #### Simplify geom
                #logger.debug("Simplification tolerance: %.10f" %tolerance)
                if geom is not None:
                    geom2  = geom.Simplify(tolerance)
            
            
        ds = None

    return geom2,xs,ys 

    
def findVertices(xoff, yoff, xsize, ysize, band, nd):
    line = band.ReadAsArray(xoff,yoff,xsize,ysize,xsize,ysize)
    if line is not None:
        nz = numpy.flatnonzero(line != nd)
    
        nzbool = nz.size > 0
        nzmin = nz[0] if nz.size > 0 else 0
        nzmax = nz[-1] if nz.size > 0 else 0
        
        return(nzbool, nzmin, nzmax)
    
    else:
        return (False,0,0)
    
    
def pl2xy(gtf,band,p,l):
    
    cols = band.XSize
    rows = band.YSize
    
    cellSizeX = gtf[1]
    cellSizeY = -1 * gtf[5]
  
    minx = gtf[0]
    maxy = gtf[3]
    
    # calc locations of pixels
    x = cellSizeX * p + minx
    y = maxy - cellSizeY * l - cellSizeY * 0.5
    
    return x,y
 
 
def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step


def buffernum(num,buf):
    sNum = str(num)
    while len(sNum)<buf:
        sNum = "0%s" %sNum
    return sNum
   
    
def deleteTempFiles(names):
    logger.info('Deleting Temp Files')
    for name in names:
        if name is not None:
            deleteList = glob.glob(os.path.splitext(name)[0]+'.*')
            for f in deleteList:
                try:
                    os.remove(f)
                    loger.info('Deleted '+os.path.basename(f))
                except:
                    logger.info('Could not remove '+os.path.basename(f))
   
                    
def copyall(srcfile,dstdir):
    for fpi in glob.glob("%s.*" %os.path.splitext(srcfile)[0]):
        fpo = os.path.join(dstdir,os.path.basename(fpi))
        shutil.copy2(fpi,fpo)
    
    
def ExecCmd(cmd):
    logger.info(cmd)
    p = Popen(cmd,shell=True,stderr=PIPE,stdout=PIPE)
    (so,se) = p.communicate()
    rc = p.wait()
    logger.info(rc)
    logger.info(se)
    logger.info(so)
    
     


def ExecCmd_mp(job):
    job_name, cmd = job
    logger.info('Running job: {0}'.format(job_name))
    logger.debug('Cmd: {0}'.format(cmd))
    
    p = Popen(cmd,shell=True,stderr=PIPE,
              stdout=PIPE,preexec_fn=os.setsid)
    try:
        (so,se) = p.communicate()
    except KeyboardInterrupt:
        os.killpg(p.pid, signal.SIGTERM)
    
    else:
        logger.debug(so)
        logger.debug(se)

def getGEMetadataAsXml(metafile):
	if os.path.isfile(metafile):
		try:
			metaf = open(metafile, "r")
		except IOError, err:
			LogMsg("Could not open metadata file %s because %s" % (metafile, err))
			raise
	else:
		LogMsg("Metadata file %s not found" % metafile)
		return None

	# Patterns to extract tag/value pairs and BEGIN/END group tags
	gepat1 = re.compile(r'(?P<tag>\w+) = "?(?P<data>.*?)"?;', re.I)
	gepat2 = re.compile(r"(?P<tag>\w+) = ", re.I)

	# These tags use the following tag/value as an attribute of the group rather than
	# a standalone node
	group_tags = {"aoiGeoCoordinate":"coordinateNumber",
				  "aoiMapCoordinate":"coordinateNumber",
				  "bandSpecificInformation":"bandNumber"}

	# Start processing
	root = ET.Element("root")
	parent = None
	current = root
	node_stack = []
	mlstr = False  # multi-line string flag

	for line in metaf:
		# mlstr will be true when working on a multi-line string
		if mlstr:
			if not line.strip() == ");":
				data += line.strip()
			else:
				data += line.strip()
				child = ET.SubElement(current, tag)
				child.text = data
				mlstr = False

		# Handle tag/value pairs and groups
		mat1 = gepat1.search(line)
		if mat1:
			tag = mat1.group("tag").strip()
			data = mat1.group("data").strip()

			if tag == "BEGIN_GROUP":
				if data is None or data == "":
					child = ET.SubElement(current, "group")
				else:
					child = ET.SubElement(current, data)
				if parent:
					node_stack.append(parent)
				parent = current
				current = child
			elif tag == "END_GROUP":
				current = parent if parent else root
				parent = node_stack.pop() if node_stack else None
			else:
				if current.tag in group_tags and tag == group_tags[current.tag]:
					current.set(tag, data)
				else:
					child = ET.SubElement(current, tag)
					child.text = data
		else:
			mat2 = gepat2.search(line)
			if mat2:
				tag = mat2.group("tag").strip()
				data = ""
				mlstr = True

	metaf.close()
	#print ET.ElementTree(root)
	return ET.ElementTree(root)


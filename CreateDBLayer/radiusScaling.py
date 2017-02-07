from xml.dom import minidom
import sys

scaling=float(sys.argv[1])/100.

def changer1(inc,oldelement):
    newstring=oldelement.nodeValue
    splittedstring=newstring.split(' ')
    splittedstring[1]=str(float(splittedstring[1])*(1+inc))
    #print newstring
    newstring=''
    cnt=0
    for tmps in splittedstring:
       if cnt==0:
          newstring=newstring+tmps
       else:
          newstring=newstring+" "+tmps
       cnt=cnt+1
       
    #print newstring
    return newstring

xmldoc = minidom.parse('HPD.xml')
itemlist = xmldoc.getElementsByTagName('paramVector') 
#print len(itemlist)
#print itemlist[0].attributes['name'].value
for s in itemlist :
    issim=s.getAttribute('name')
    #print ciccio
#    changer1(1.0,s.firstChild)
    if 'sim' not in issim:
       s.firstChild.nodeValue=changer1(scaling,s.firstChild)
    #help(s)
    #for sc in s.childNodes
    #print sc.attributes['name'].value
    
    #print s.toxml()

#final=xmldoc.toxml()
file_handle = open("HPDscaling_%i_perc.xml" % int(scaling*100),"wb")
xmldoc.writexml(file_handle)
file_handle.close()

#print final

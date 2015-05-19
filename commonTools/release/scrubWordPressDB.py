#! /usr/bin/env python

import xml.etree.ElementTree as xml
import csv

namespaces = {"excerpt" : "http://wordpress.org/export/1.2/excerpt/",
              "content" : "http://purl.org/rss/1.0/modules/content/",
              "wfw"     : "http://wellformedweb.org/CommentAPI/",
              "dc"      : "http://purl.org/dc/elements/1.1/",
              "wp"      : "http://wordpress.org/export/1.2/"}

wpPath      = "thetrilinosproject.wordpress.2015-04-30.xml"
packagePath = "../../PackagesList.cmake"
outPath     = "Descriptions.csv"

packages = []
data = open(packagePath, 'r').readlines()
start = data.index("TRIBITS_REPOSITORY_DEFINE_PACKAGES(\n")
stop  = data.index("  )\n")
for line in data[start+1:stop]:
    columns = line.split()
    packages.append(columns[0])

packages_low = [pkg.lower() for pkg in packages]

tree = xml.parse(wpPath)
root = tree.getroot()

descriptions = []
for item in root.iter("item"):
    title = item.find("title").text
    if title:
        title_low = title.lower()
        try:
            index = packages_low.index(title_low)
            description = item.find("content:encoded", namespaces).text
            if description:
                packageID = packages[index] + "-Devel"
                description = description.encode('ascii', 'ignore')
                packageData = {"package_id"  : packageID,
                               "description" : description}
                descriptions.append(packageData)
        except ValueError:
            pass

fieldnames = ["package_id","description"]
with open(outPath, 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for packageData in descriptions:
        writer.writerow(packageData)

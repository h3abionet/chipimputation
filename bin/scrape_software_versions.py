#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'h3ahipimputation': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_vcftools.txt', r"vcftools v(\S+)"],
    'plink2': ['v_plink2.txt', r"plink2, version (\S+)"],
}
results = OrderedDict()
results['h3achipimputation'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['vcftools'] = '<span style="color:#999999;\">N/A</span>'
results['plink2'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'h3achipimputation-software-versions'
section_name: 'h3achipimputation Software Versions'
section_href: 'https://github.com/h3abionet/chipimputation'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")

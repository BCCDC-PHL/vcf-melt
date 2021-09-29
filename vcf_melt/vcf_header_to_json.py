#!/usr/bin/env python

import argparse
import json
import re


def parse_untagged_record(record):
    parsed_record = {}
    [k, v] = record.split('=', 1)
    k = k.strip()
    v = v.strip('"').strip()
    parsed_record[k] = v
    return parsed_record


def parse_tag(tag):
    parsed_tag = {}
    tag = tag.strip('<').strip('>')
    quoted_tags = re.findall('\w+=\"[^\"]+\"', tag)
    unquoted_tags = re.findall('\w+=[^\",]+', tag)
    unquoted_kv_pairs = [e.split('=') for e in unquoted_tags]
    quoted_kv_pairs = [e.split('=') for e in quoted_tags]
    for pair in unquoted_kv_pairs + quoted_kv_pairs:
        pair[1]= pair[1].strip('"')
        parsed_tag[pair[0]] = pair[1]

    return parsed_tag


def parse_tagged_record(record):
    parsed_record = {}
    [top_level_key, tag] = record.split('=', 1)
    parsed_tag = parse_tag(tag)
    parsed_record[top_level_key] = parsed_tag
    return parsed_record



def main(args):
    output = {}
    with open(args.vcf, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('##'):
                line = line.strip('#')
                parsed_record = {}
                # print(line)
                if re.match("\w+=<[^<>]+>", line):
                    parsed_record = parse_tagged_record(line)
                elif re.match("\w+=*", line):
                    parsed_record = parse_untagged_record(line)
                #print(json.dumps(parsed_record))
                #continue
                for k, v in parsed_record.items():
                    if k in output:
                        if type(output[k]) == type([]):
                            output[k].append(v)
                        else:
                            output[k] = [output[k], v]
                    else:
                        output[k] = v

    print(json.dumps(output, indent=2))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf')
    args = parser.parse_args()
    main(args)

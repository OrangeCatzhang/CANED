#!/usr/bin/env python
import os
import requests


def get_ec_number(pdb_id):
    # GraphQL 查询
    query = """
    query GetEnzymeECByPDBID($pdbid: String!) {
      polymer_entity(entry_id: $pdbid, entity_id: "1") {
        rcsb_polymer_entity {
          pdbx_ec
        }
      }
    }
    """
    variables = {
        "pdbid": pdb_id
    }
    headers = {
        "Content-Type": "application/json"
    }
    url = "https://data.rcsb.org/graphql"

    try:
        response = requests.post(url, json={"query": query, "variables": variables}, headers=headers)
        response.raise_for_status()  # 将引发 HTTPError，如果状态不是 200
        data = response.json()
        if "data" in data and data["data"]["polymer_entity"]:
            enzyme_ec = data["data"]["polymer_entity"]["rcsb_polymer_entity"]["pdbx_ec"]
            print(f"PDB ID: {pdb_id}, EC number: {enzyme_ec}")
            return enzyme_ec
        else:
            print(f"No data found for PDB ID {pdb_id}")
    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except requests.exceptions.ConnectionError as conn_err:
        print(f"Connection error occurred: {conn_err}")
    except requests.exceptions.Timeout as timeout_err:
        print(f"Timeout error occurred: {timeout_err}")
    except requests.exceptions.RequestException as req_err:
        print(f"An error occurred: {req_err}")
    return None


def create_index(pdb_directory, index_file="index.txt", last_processed_file="last_processed.txt"):
    last_processed_pdb = None
    if os.path.exists(last_processed_file):
        with open(last_processed_file, "r") as f:
            last_processed_pdb = f.read().strip()

    # 使用绝对路径保存索引文件和进度文件
    index_file = os.path.abspath(index_file)
    last_processed_file = os.path.abspath(last_processed_file)

    # 检查索引文件是否已存在，并读取已有的 PDB ID
    existing_pdb_ids = []
    if os.path.exists(index_file):
        with open(index_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) > 0:
                    existing_pdb_ids.append(parts[0])

    # 打开索引文件以追加模式写入
    with open(index_file, "a") as index_f, open(last_processed_file, "w") as last_f:
        for pdb_id in sorted(os.listdir(pdb_directory)):
            # 如果之前已经处理过此 PDB ID 或者已被索引，跳过
            if (last_processed_pdb and pdb_id <= last_processed_pdb) or pdb_id in existing_pdb_ids:
                continue

            pdb_path = os.path.join(pdb_directory, pdb_id)
            if os.path.isdir(pdb_path):
                print(f"Processing {pdb_id} in {pdb_path}")

                # 查询 EC 编号
                ec_number = get_ec_number(pdb_id)
                if ec_number:
                    index_f.write(f"{pdb_id}\t{ec_number}\t{pdb_path}\n")
                else:
                    print(f"EC number not found for PDB ID {pdb_id}, skipping.")
                    with open("failed_pdb.txt", "a") as failed_f:
                        failed_f.write(f"{pdb_id}\n")
                last_f.seek(0)  # 移动到文件开头
                last_f.truncate()  # 清空文件内容
                last_f.write(pdb_id)  # 写入当前 PDB ID

    print(f"Index file saved to: {index_file}")

def search_by_pdb(index_file, pdb_id):
    """
    根据PDB名称搜索索引文件。
    :param index_file: 索引文件路径
    :param pdb_id: PDB名称
    :return: 匹配的目录路径
    """
    with open(index_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts[0] == pdb_id:
                return parts[2]
    return None


import re

def search_by_ec(index_file, ec_number):
    """
    根据EC编号搜索索引文件, 支持精确匹配和模糊匹配。
    :param index_file: 索引文件路径
    :param ec_number: EC编号
    :return: 匹配的目录路径列表
    """
    result = []
    # 检查EC编号是否为精确匹配（四个部分）
    is_exact_match = re.match(r"^\d+\.\d+\.\d+\.\d+$", ec_number)
    
    # 构建正则表达式
    if is_exact_match:
        # 精确匹配
        pattern = re.compile(rf"^{re.escape(ec_number)}$")
    else:
        # 模糊匹配
        pattern = re.compile(rf"{re.escape(ec_number)}")

    with open(index_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 3:
                continue  
            pdbid, current_ec, path = parts
            if pattern.search(current_ec):  
                result.append(path)
    return result



import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PDB数据库索引和搜索工具")
    parser.add_argument("-p","--pdb_directory", type=str, help="PDB文件目录路径")
    parser.add_argument("-i","--index_file", type=str, default="index.txt", help="索引文件路径")
    parser.add_argument("--last_processed_file", type=str, default="last_processed.txt", help="记录最后处理的文件路径")
    parser.add_argument("-sp","--search_pdb", type=str, help="根据PDB ID搜索索引文件")
    parser.add_argument("-se","--search_ec", type=str, help="根据EC编号搜索索引文件")

    args = parser.parse_args()

    if args.pdb_directory:
        create_index(args.pdb_directory, args.index_file, args.last_processed_file)

    if args.search_pdb:
        pdb_path = search_by_pdb(args.index_file, args.search_pdb)
        if pdb_path:
            print(f"PDB ID {args.search_pdb} found at: {pdb_path}")
        else:
            print(f"PDB ID {args.search_pdb} not found.")

    if args.search_ec:
        ec_paths = search_by_ec(args.index_file, args.search_ec)
        if ec_paths:
            print(f"EC Number {args.search_ec} found at:")
            for path in ec_paths:
                print(path)
        else:
            print(f"EC Number {args.search_ec} not found.")
#!/usr/bin/env python3

import argparse
import pyrosetta
import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts
import pandas as pd


def relax_with_constraints(pdb_path, cst_file, output_path, cutoff=0):
    """
    Use Rosetta FastRelax protocol to optimize a PDB structure with constraints.
    
    Parameters:
        pdb_path (str): Path to input PDB file
        cst_file (str): Path to constraint file (.cst)
        output_path (str): Path to output optimized PDB file
        cutoff (float): Energy cutoff for cst_score filter, default is 0 (no filtering if 0)
    
    Returns:
        pose: Optimized structure
        df_scores: DataFrame of scores
    """
    # Load PDB file
    pose = pyrosetta.pose_from_pdb(pdb_path)
    print(f"Loaded PDB: {pdb_path}")
    
    # 根据 cutoff 设置 confidence 值：cutoff=0 时为 1（不筛选），否则为 0（启用筛选）
    confidence = "0" if cutoff == 0 else "1"
    
    # Define RosettaScripts XML with dynamic cutoff and confidence
    xml = f"""
    <ROSETTASCRIPTS>  
      <SCOREFXNS>
          <ScoreFunction name="sfxn_design" weights="beta">
              <Reweight scoretype="atom_pair_constraint" weight="1.0"/>
              <Reweight scoretype="dihedral_constraint" weight="1.0"/>
              <Reweight scoretype="angle_constraint" weight="1.0"/>
              <Reweight scoretype="coordinate_constraint" weight="1.0"/>
          </ScoreFunction>
          <ScoreFunction name="sfxn" weights="beta"/>
      </SCOREFXNS>
      <RESIDUE_SELECTORS>
          <Chain name="chainA" chains="A"/>
          <Not name="chainB" selector="chainA"/>
          <Or name="everything" selectors="chainA,chainB"/>
      </RESIDUE_SELECTORS>
      <TASKOPERATIONS>
          <LimitAromaChi2 name="limitchi2" chi2max="110" chi2min="70" include_trp="True" />
          <ExtraRotamersGeneric name="ex1_ex2aro" ex1="1" ex2aro="1" ex2="0"/>
          <IncludeCurrent name="ic"/>
          <OperateOnResidueSubset name="to_relax" selector="everything"> 
              <RestrictToRepackingRLT/>
          </OperateOnResidueSubset>
      </TASKOPERATIONS>
      <MOVERS>    
          <AddOrRemoveMatchCsts name="add_enz_csts" cstfile="{cst_file}" cst_instruction="add_new"/>
          <FastRelax name="fastRelax" scorefxn="sfxn_design" repeats="3" task_operations="ex1_ex2aro,ic,limitchi2,to_relax" relaxscript="MonomerRelax2019">
              <MoveMap name="MM" bb="1" chi="1" jump="1"/>
          </FastRelax>
      </MOVERS>
      <FILTERS>
          <EnzScore name="cst_score" scorefxn="sfxn_design" confidence="{confidence}" whole_pose="1" score_type="cstE" energy_cutoff="{cutoff}"/>
      </FILTERS>
      <PROTOCOLS>
          <Add mover="add_enz_csts"/>
          <Add mover="fastRelax"/>
          <Add filter="cst_score"/>
      </PROTOCOLS>
    </ROSETTASCRIPTS>
    """
    
    # Create and run RosettaScripts task
    try:
        task_relax = rosetta_scripts.SingleoutputRosettaScriptsTask(xml)
        task_relax.setup()  # Check XML syntax
        print("Starting FastRelax with constraints...")
        packed_pose = task_relax(pose)
        print("Relax completed.")

        # Check if cst_score passes the filter (only if cutoff != 0)
        
        if packed_pose.scores is not None:
            # Save optimized structure
            output_pdb = f"{output_path}.pdb"
            packed_pose.pose.dump_pdb(output_pdb)
            print(f"Saved optimized structure to: {output_pdb}")
            # Save scores to CSV file
            df_scores = pd.DataFrame.from_records([packed_pose.scores])
            output_scores = f"{output_path}_scores.sc"
            df_scores.to_csv(output_scores, sep=' ', index=False, lineterminator='\n')
            print(f"Saved scores to: {output_scores}")
              
        return packed_pose.pose, df_scores
    except Exception as e:
        print(f"Error during FastRelax: {e}")
        return None, None

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Perform FastRelax with constraints on a PDB file.")
    parser.add_argument("-i", "--pdb", type=str, required=True, help="Path to input PDB file")
    parser.add_argument("-c", "--cst", type=str, required=True, help="Path to constraint file (.cst)")
    parser.add_argument("-o", "--output", type=str, default="relaxed", help="Path to output relaxed PDB file")
    parser.add_argument("-pa", "--params", type=str, required=True, help="Path to params file for extra_res_fa")
    parser.add_argument("-cut", "--cutoff", type=float, default=0, help="Energy cutoff for cst_score filter (0 = no filtering)")
    
    args = parser.parse_args()
    
    pyrosetta.init(extra_options=f"-run:preserve_header -beta -extra_res_fa {args.params} -ignore_unrecognized_res false -ignore_zero_occupancy false")
    
    # Run Relax
    pose, df_scores = relax_with_constraints(args.pdb, args.cst, args.output, args.cutoff)

if __name__ == "__main__":
    main()
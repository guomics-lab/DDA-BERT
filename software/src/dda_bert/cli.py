from __future__ import annotations

import multiprocessing

from . import common_config

config_data = common_config.load_base_config()

import typer
from rich.console import Console
from pathlib import Path
from typing import Literal, Optional, Annotated
from importlib.metadata import version

from .obj.Entity import InputParam
from .utils import file_utils

MassFileFormat = Literal["raw", "wiff", "d"]
SearchEngines = Literal["sage", "sage,fp", "sage,fp,ap"]

app = typer.Typer(
    name="DDA-BERT",
    help="An end-to-end rescoring tool tailored for data-dependent acquisition (DDA) proteomics.",
    no_args_is_help=True,
    add_completion=False,
)

console = Console()


def _version_callback(value: bool):
    if value:
        __version__ = version("dda-bert")
        print(f"Version: {__version__}")
        raise typer.Exit()


def _get_mzml_list(mzml_paths: Optional[str], txt_path: Optional[Path]) -> list[str]:
    if mzml_paths:
        return mzml_paths.split(',')
    elif txt_path:
        return file_utils.read_txt(txt_path)
    else:
        raise typer.BadParameter("One of --mzml-paths or --txt-path must be provided.")


def _prepare_gpu_devices(gpu_devices: str) -> tuple[list[int], int]:
    from .utils.gpu_utils import get_usage_device, get_top_free_max_device
    if gpu_devices == 'auto':
        max_member_use_rate = config_data.get('gpu', {}).get('max_member_use_rate', 0.5)
        usage_device_list, _ = get_usage_device(max_member_use_rate)
        input_device = get_top_free_max_device(usage_device_list)
        console.print(f'Use device list: {usage_device_list}, max free member device is: {input_device}')
        return usage_device_list, input_device
    else:
        try:
            usage_device_list = [int(nn) for nn in gpu_devices.split(',')]
            input_device = get_top_free_max_device(usage_device_list)
            console.print(f'Use device list: {usage_device_list}, max free member device is: {input_device}')
            return usage_device_list, input_device
        except ValueError:
            raise typer.BadParameter(
                "Invalid format for --gpu-devices. Must be comma-separated integers (e.g., '0,1').")


def _prepare_input_param(params: dict, only_run_predict: bool) -> InputParam:
    input_param = InputParam()
    input_param.env = 'linux'
    input_param.only_run_predict = only_run_predict

    input_param.thread_num = params['thread_num']
    input_param.mass_file_format = params['mass_file_format']
    input_param.search_engines = params['search_engines']
    input_param.output_dir_path = params['output_path']
    input_param.model_path = params['model_path']
    input_param.mzml_file_list = params['mzml_file_list']

    input_param.gpu_devices, input_param.device = _prepare_gpu_devices(params['gpu_devices'])

    if not only_run_predict:
        input_param.fasta_path = params['fasta']
        input_param.sage_config_file = params['sage_config_file']
        input_param.fp_workflow_file = params['fp_workflow_file']
        input_param.ap_config_file = params['ap_config_file']
    else:
        input_param.mzml_ipc_file_dir = params.get('mzml_ipc_file_dir')
        input_param.sage_file_dir = params.get('sage_file_dir')
        input_param.fp_file_dir = params.get('fp_file_dir')
        input_param.ap_file_dir = params.get('ap_file_dir')
        input_param.fasta_path = params.get('fasta')

    return input_param


@app.callback()
def main_version(
        version: bool = typer.Option(
            None,
            "--version",
            callback=_version_callback,
            is_eager=True,
            help="Show the version and exit.",
        ),
):
    pass


@app.command(
    help="Full assessment workflow: First perform identification using the specified search engine(s), then conduct quality scoring.")
def assess(
        mzml_paths: Annotated[
            Optional[str], typer.Option("--mzml-paths", help="mzML file absolute path list, split by ',' ")] = None,
        txt_path: Annotated[Optional[str], typer.Option("--txt-path",
                                                        help="mzML file absolute path list txt file path. Each mzml in one line")] = None,
        output_path: Annotated[str, typer.Option("--output", "-o",
                                                 help="Path where search and quant results will be written. (Required)")] = ...,
        fasta: Annotated[str, typer.Option("--fasta", "-f", help="Path to FASTA database. (Required)")] = ...,

        thread_num: Annotated[int, typer.Option("--thread-num", help="Run thread number.")] = 20,
        model_path: Annotated[str, typer.Option("--model", help="Predict model path.")] = Path(
            config_data['model_file']),
        gpu_devices: Annotated[
            str, typer.Option("--gpu-devices", help="GPU devices index list. Split by ',' or 'auto'")] = 'auto',
        mass_file_format: Annotated[str, typer.Option("--mass-format",
                                                      help="The type of the raw mass spectrometry file: raw, wiff, d")] = "raw",
        search_engines: Annotated[
            str, typer.Option("--engines", help="Search engine combination: sage, fp, ap")] = "sage,fp,ap",

        sage_config_file: Annotated[Optional[str], typer.Option("--sage-config-file",
                                                                help="Sage run config, if not set will run use default config")] = None,
        fp_workflow_file: Annotated[Optional[str], typer.Option("--fp-workflow-file",
                                                                help="Fragpipe run workflow, if not set will run use default config")] = None,
        ap_config_file: Annotated[Optional[str], typer.Option("--ap-config-file",
                                                              help="Alphapept run config, if not set will run use default config")] = None, ):
    console.print(
        'Full assessment workflow: First perform identification using the specified search engine(s), then conduct quality scoring.')
    mzml_list = _get_mzml_list(mzml_paths, txt_path)
    params = locals().copy()
    params['mzml_file_list'] = mzml_list
    input_param = _prepare_input_param(params, only_run_predict=False)

    from .main_process_handler import MainProcessHandler
    MainProcessHandler(input_param)._processing()


@app.command(
    help="Score only, skip identification (suitable for data with existing identification results, significantly faster)")
def score(
        mzml_paths: Annotated[
            Optional[str], typer.Option("--mzml-paths", help="mzML file absolute path list, split by ',' ")] = None,
        txt_path: Annotated[Optional[str], typer.Option("--txt-path",
                                                        help="mzML file absolute path list txt file path. Each mzml in one line")] = None,
        output_path: Annotated[str, typer.Option("--output", "-o",
                                                 help="Path where search and quant results will be written. (Required)")] = ...,
        fasta: Annotated[
            Optional[str], typer.Option("--fasta", "-f", help="Path to FASTA database for context. (Required)")] = ...,

        mzml_ipc_file_dir: Annotated[Optional[str], typer.Option("--mzml-ipc-file-dir",
                                                                 help="The mzml info file path, file name must be like: xxx.ipc.")] = None,
        sage_file_dir: Annotated[Optional[str], typer.Option("--sage-file-dir",
                                                             help="The sage result tsv file dir path, the file name must be: xxx.sage.tsv. Must be get if search engine contains sage.")] = None,
        fp_file_dir: Annotated[Optional[str], typer.Option("--fp-file-dir",
                                                           help="The fragpipe result pin file path. Must be get if search engine contains fp.")] = None,
        ap_file_dir: Annotated[Optional[str], typer.Option("--ap-file-dir",
                                                           help="The alphapept result hdf file path. Must be get if search engine contains ap.")] = None,

        thread_num: Annotated[int, typer.Option("--thread-num", help="Run thread number.")] = 20,
        model_path: Annotated[str, typer.Option("--model", help="Predict model path.")] = Path(
            config_data['model_file']),
        gpu_devices: Annotated[
            str, typer.Option("--gpu-devices", help="GPU devices index list. Split by ',' or 'auto'")] = 'auto',
        mass_file_format: Annotated[str, typer.Option("--mass-format",
                                                      help="The type of the raw mass spectrometry file: raw, wiff, d")] = "raw",
        search_engines: Annotated[
            str, typer.Option("--engines",
                              help="Search engine combination. Can sage or sage,fp or sage,fp,ap. Default is sage,fp,ap")] = "sage,fp,ap",
):
    console.print(
        'Score only, skip identification (suitable for data with existing identification results, significantly faster)')

    mzml_list = _get_mzml_list(mzml_paths, txt_path)

    engines = {e.strip() for e in search_engines.split(",")}
    if "sage" in engines and not sage_file_dir:
        raise typer.BadParameter("--engines contains sage, but --sage-file-dir is missing.")
    if "fp" in engines and not fp_file_dir:
        raise typer.BadParameter("--engines contains fp, but --fp-file-dir is missing.")
    if "ap" in engines and not ap_file_dir:
        raise typer.BadParameter("--engines contains ap, but --ap-file-dir is missing.")

    params = locals().copy()
    params['mzml_file_list'] = mzml_list
    input_param = _prepare_input_param(params, only_run_predict=True)

    from .main_process_handler import MainProcessHandler
    MainProcessHandler(input_param)._processing()


def main():
    app()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    multiprocessing.set_start_method('spawn', force=True)
    main()

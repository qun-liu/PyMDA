#!/usr/bin/env python3
"""Launch the local-only PyMDA web interface."""
import argparse
import uvicorn

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port',type=int,default=8000)
    args = parser.parse_args()
    uvicorn.run('main:app',host='127.0.0.1',port=args.port,proxy_headers=False)
